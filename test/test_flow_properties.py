#!/usr/bin/env python3
"""
src_helper/flow_properties（軌道に沿った流れの量）のテスト。

ソルバーは走らせないツールなので、確かめるのは式と読み書き:

  1. Mach 数は、ソルバーが Mach の空力表で使う式（satellite.get_mach_number）と一致する。
     比熱比と分子量は制御ファイル（config.yml の atmosphere 節）から読む
  2. 設定の優先順位: コマンドライン > 制御ファイル > 気体の既定値
  3. Reynolds 数の粘性係数は Sutherland 式（T0 で mu0）、圧力は rho R T
  4. 熱流束は Sutton-Graves の式の値そのもの、Tauber の式は W/m2（v1 は kW/m2 のまま
     W/m2 と書いていた）、総加熱量は熱流束の時間積分
  5. 書き出したファイルは Tecplot の読み取り（animate_trajectory/tecplot_reader）で読み直せる。
     見出しの変数名は値の列と数が合う（v1 はカンマが抜け、総加熱量の見出しが無かった）
  6. 風があれば対気速度を使い、Kn はソルバーの列をそのまま写す
"""

import argparse
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml

from context import ROOT_DIR

HELPER_DIR = os.path.join(ROOT_DIR, 'src_helper')
for name in ('flow_properties', 'animate_trajectory'):
    directory = os.path.join(HELPER_DIR, name)
    if directory not in sys.path:
        sys.path.insert(0, directory)

import flow_properties as flow_properties  # noqa: E402
import tecplot_reader as tecplot_reader    # noqa: E402
import satellite.satellite as satellite    # noqa: E402

SCRIPT = os.path.join(HELPER_DIR, 'flow_properties', 'flow_properties.py')


def write_tacode_output(path, wind=False):
    """Tacode の tecplot.dat と同じ列を持つ小さなファイル。"""
    time = np.arange(0.0, 5.0)
    speed = np.array([6000.0, 5900.0, 5700.0, 5300.0, 4600.0])
    density = np.array([1.e-6, 3.e-6, 1.e-5, 3.e-5, 1.e-4])
    temperature = np.array([150.0, 155.0, 160.0, 170.0, 185.0])
    knudsen = np.array([1.e-1, 3.e-2, 1.e-2, 3.e-3, 1.e-3])
    variables = ['Time[s]', 'X[km]', 'Y[km]', 'Z[km]', 'Long[deg.]', 'Lati[deg.]', 'Alti[km]', 'Upl[m/s]',
                 'Vpl[m/s]', 'Wpl[m/s]', 'VelplAbs[m/s]', 'Dens[kg/m3]', 'Temp[K]', 'Kn']
    if wind:
        variables += ['WindE[m/s]', 'WindN[m/s]', 'WindU[m/s]', 'VelairAbs[m/s]']
    with open(path, 'w') as f:
        f.write('# Tecplot data: Tacode\n')
        f.write('Variables = ' + ','.join(variables) + '\n')
        f.write('zone t=time i= {} f=point\n'.format(len(time)))
        for n in range(len(time)):
            row = [time[n], 0.0, 0.0, 0.0, 0.0, 0.0, 80.0 - 10.0*n, speed[n], 0.0, 0.0, speed[n],
                   density[n], temperature[n], knudsen[n]]
            if wind:
                row += [30.0, 0.0, 0.0, speed[n] - 30.0]
            f.write(' '.join(repr(float(v)) for v in row) + '\n')
    return {'time': time, 'speed': speed, 'density': density, 'temperature': temperature, 'knudsen': knudsen}


def namespace(**change):
    value = {'control': None, 'gas': 'co2', 'specific_heat_ratio': None, 'molecular_weight': None,
             'length': 2.0, 'sutherland': None, 'nose_radius': 0.5, 'heating': 'sutton-graves',
             'heat_flux_law': None}
    value.update(change)
    return argparse.Namespace(**value)


class TemporaryDirectoryMixin:

    def setUp(self):
        self.directory_object = tempfile.TemporaryDirectory()
        self.directory = self.directory_object.name

    def tearDown(self):
        self.directory_object.cleanup()


class TestTheSettings(TemporaryDirectoryMixin, unittest.TestCase):

    def write_control(self, atmosphere, satellite_section=None):
        path = os.path.join(self.directory, 'config.yml')
        with open(path, 'w') as f:
            yaml.safe_dump({'atmosphere': atmosphere, 'satellite': satellite_section or {}}, f)
        return path

    def test_the_control_file_gives_the_gas_and_the_length(self):
        path = self.write_control({'specific_heat_ratio': 1.335, 'molecular_weight': 43.49},
                                  {'characteristic_length': 1.325})
        value, source = flow_properties.resolve_setting(namespace(control=path, length=None))
        self.assertEqual((value['specific_heat_ratio'], value['molecular_weight'], value['length']),
                         (1.335, 43.49, 1.325))
        self.assertEqual(source['molecular_weight'], path)

    def test_the_command_line_wins_over_the_control_file(self):
        path = self.write_control({'specific_heat_ratio': 1.335, 'molecular_weight': 43.49})
        value, source = flow_properties.resolve_setting(namespace(control=path, specific_heat_ratio=1.3))
        self.assertEqual(value['specific_heat_ratio'], 1.3)
        self.assertEqual(source['specific_heat_ratio'], 'given')
        self.assertEqual(value['molecular_weight'], 43.49)

    def test_without_them_the_gas_gives_the_default(self):
        value, _ = flow_properties.resolve_setting(namespace(gas='air'))
        self.assertEqual((value['specific_heat_ratio'], value['molecular_weight']), (1.4, 28.9644))

    def test_a_missing_length_stops(self):
        with self.assertRaises(SystemExit):
            flow_properties.resolve_setting(namespace(length=None))

    def test_a_missing_nose_radius_is_reported_and_the_placeholder_used(self):
        import contextlib
        import io
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            value, source = flow_properties.resolve_setting(namespace(nose_radius=None))
        self.assertEqual(value['nose_radius'], flow_properties.NOSE_RADIUS_DEFAULT)
        self.assertEqual(source['nose_radius'], 'not given')
        self.assertIn('nose radius is not given', buffer.getvalue())

    def test_a_given_nose_radius_says_nothing(self):
        import contextlib
        import io
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            value, _ = flow_properties.resolve_setting(namespace(nose_radius=0.22))
        self.assertEqual(value['nose_radius'], 0.22)
        self.assertNotIn('nose radius', buffer.getvalue())

    def test_custom_heating_without_its_law_stops(self):
        with self.assertRaises(SystemExit):
            flow_properties.resolve_setting(namespace(heating='custom'))


class TestTheFormulas(TemporaryDirectoryMixin, unittest.TestCase):

    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.directory, 'tecplot.dat')
        self.given = write_tacode_output(self.path)
        self.data = tecplot_reader.read_tecplot(self.path)

    def compute(self, **change):
        value, _ = flow_properties.resolve_setting(namespace(**change))
        return flow_properties.compute(self.data, value, value['nose_radius']), value

    def test_the_mach_number_is_the_solvers(self):
        result, value = self.compute(specific_heat_ratio=1.335, molecular_weight=43.49)
        factor = 1.335*satellite.CONSTANT_GAS_UNIVERSAL/(43.49*1.e-3)
        expected = [satellite.get_mach_number(v, t, {satellite.KEY_SOUND_SPEED: factor})
                    for v, t in zip(self.given['speed'], self.given['temperature'])]
        np.testing.assert_allclose(result['Mach'], expected, rtol=1.e-15)

    def test_the_gas_constant_is_the_solvers(self):
        self.assertEqual(flow_properties.CONSTANT_GAS_UNIVERSAL, satellite.CONSTANT_GAS_UNIVERSAL)

    def test_the_pressure_is_rho_r_t(self):
        result, value = self.compute()
        constant = flow_properties.CONSTANT_GAS_UNIVERSAL/(value['molecular_weight']*1.e-3)
        np.testing.assert_allclose(result['Pres[Pa]'], self.given['density']*constant*self.given['temperature'],
                                   rtol=1.e-15)

    def test_the_viscosity_is_mu0_at_t0(self):
        mu0, t0, constant = flow_properties.GAS['co2']['sutherland']
        self.assertAlmostEqual(flow_properties.viscosity_sutherland(t0, mu0, t0, constant), mu0, places=20)
        mu0, t0, constant = flow_properties.GAS['air']['sutherland']
        # 空気の 300 K で 1.846e-5 Pa s（White, Viscous Fluid Flow 表 1-2 の Sutherland 式と同じ値）
        self.assertAlmostEqual(flow_properties.viscosity_sutherland(300.0, mu0, t0, constant), 1.846e-5, delta=0.005e-5)

    def test_the_reynolds_number(self):
        result, value = self.compute(length=2.0)
        viscosity = flow_properties.viscosity_sutherland(self.given['temperature'], *value['sutherland'])
        np.testing.assert_allclose(result['Re'], self.given['density']*self.given['speed']*2.0/viscosity,
                                   rtol=1.e-15)

    def test_sutton_graves(self):
        result, _ = self.compute(nose_radius=0.5)
        expected = 1.9027e-4*np.sqrt(self.given['density']/0.5)*self.given['speed']**3
        np.testing.assert_allclose(result['HeatFlux[W/m2]'], expected, rtol=1.e-14)

    def test_tauber_is_in_watts_per_square_metre(self):
        # Tauber の火星の式と Sutton-Graves は同じ量を数 % で与える（v1 の係数 1.35e-7 だと
        # 3 桁小さい kW/m2 になる）
        tauber, _ = self.compute(heating='tauber-mars')
        sutton, _ = self.compute()
        ratio = tauber['HeatFlux[W/m2]']/sutton['HeatFlux[W/m2]']
        self.assertTrue(np.all((ratio > 0.9) & (ratio < 1.25)), ratio)

    def test_the_heat_load_is_the_integral_of_the_flux(self):
        result, _ = self.compute(heating='custom', heat_flux_law=[2.0, 0.0, 0.0], nose_radius=1.0)
        np.testing.assert_allclose(result['HeatFlux[W/m2]'], 2.0)
        np.testing.assert_allclose(result['HeatLoad[J/m2]'], 2.0*self.given['time'])

    def test_the_knudsen_number_is_the_solvers_column(self):
        result, _ = self.compute()
        np.testing.assert_array_equal(result['Kn'], self.given['knudsen'])

    def test_the_air_relative_speed_is_used_with_wind(self):
        path = os.path.join(self.directory, 'tecplot_wind.dat')
        write_tacode_output(path, wind=True)
        value, _ = flow_properties.resolve_setting(namespace())
        result = flow_properties.compute(tecplot_reader.read_tecplot(path), value, 0.5)
        np.testing.assert_allclose(result['Velair[m/s]'], self.given['speed'] - 30.0)


class TestTheScript(TemporaryDirectoryMixin, unittest.TestCase):

    def test_it_writes_a_file_the_reader_reads_back(self):
        path = os.path.join(self.directory, 'tecplot.dat')
        given = write_tacode_output(path)
        with open(os.path.join(self.directory, 'config.yml'), 'w') as f:
            yaml.safe_dump({'atmosphere': {'specific_heat_ratio': 1.335, 'molecular_weight': 43.49},
                            'satellite': {'characteristic_length': 1.325}}, f)
        completed = subprocess.run([sys.executable, SCRIPT, 'tecplot.dat', '-o', 'flow.dat', '--nose-radius', '0.5'],
                                   cwd=self.directory, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn('config.yml', completed.stdout)
        self.assertNotIn('nose radius is not given', completed.stdout)

        with open(os.path.join(self.directory, 'flow.dat')) as f:
            lines = f.readlines()
        header = [line for line in lines if line.lower().startswith('variables')][0]
        names = [name.strip().strip('"') for name in header.split('=', 1)[1].split(',')]
        self.assertEqual(names, flow_properties.VARIABLES)
        row = [line for line in lines if line[:1].isdigit()][0].split()
        self.assertEqual(len(row), len(names))

        data = tecplot_reader.read_tecplot(os.path.join(self.directory, 'flow.dat'))
        np.testing.assert_allclose(data['Time'], given['time'])
        np.testing.assert_allclose(data['Dens'], given['density'], rtol=1.e-6)
        self.assertIn('HeatLoad', data)

    def test_the_script_warns_without_a_nose_radius(self):
        write_tacode_output(os.path.join(self.directory, 'tecplot.dat'))
        completed = subprocess.run([sys.executable, SCRIPT, 'tecplot.dat', '-o', 'flow.dat', '--length', '1.0'],
                                   cwd=self.directory, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn('nose radius is not given', completed.stdout)
        with open(os.path.join(self.directory, 'flow.dat')) as f:
            self.assertIn('NOT GIVEN', f.read())

    def test_a_missing_input_stops(self):
        completed = subprocess.run([sys.executable, SCRIPT, 'no_such_file.dat'], cwd=self.directory,
                                   capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)


if __name__ == '__main__':
    unittest.main()
