#!/usr/bin/env python3
"""
ケプラー運動の保存則による積分器のテスト。

J 項・自転・空力抗力をすべて 0 にすると純粋な 2 体問題になり、
比エネルギーと角運動量が保存する。RK4 の実装が壊れるとこれが崩れる。
"""

import unittest

import numpy as np

from context import two_body_config, gravitational_parameter, quiet

import coordinate_system.coordinate_system as coordinate_system
import solver.solver as solver
from orbital.orbital import orbital


def circular_orbit_config(altitude_m=4.0e5):
    """指定高度の円軌道になる初期条件を組んだ config を返す。"""
    config = two_body_config()

    radius_equat = config['planet']['radius']
    gm = gravitational_parameter(config)

    # 赤道上、経度 0 度から出発する
    config['initial_settings']['coordinate'] = [0.0, 0.0, altitude_m * 1.0e-3]

    # 円軌道速度。config の velocity は [東, 北, 上] m/s
    radius = radius_equat + altitude_m
    speed = np.sqrt(gm / radius)
    config['initial_settings']['velocity'] = [speed, 0.0, 0.0]
    config['initial_settings']['density_factor'] = [0.0]

    return config


def integrate(config, time_max, timestep):
    """ソルバーを回して軌道と速度の履歴を返す。"""
    config['computational_setup']['time_elapsed_maximum'] = time_max
    config['time_integration']['timestep_constant'] = timestep

    with quiet():
        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(
                config, iteration, time_elapsed,
                coordinate_dict, velocity_dict, trajectory_dict,
                atmosphere_dict={}, aerodynamic_dict={})

    position = np.array(coordinate_dict['cartesian'])
    velocity = np.array(velocity_dict['cartesian'])
    return iteration, position, velocity


class TestTwoBodyConservation(unittest.TestCase):
    """2 体問題における保存量。"""

    @classmethod
    def setUpClass(cls):
        config = circular_orbit_config()
        cls.gm = gravitational_parameter(config)
        # 低軌道の周期は約 5,500 s。1 周期強を 2 s 刻みで積分する
        cls.iteration, cls.position, cls.velocity = integrate(config, 6000.0, 2.0)

    def test_specific_energy_is_conserved(self):
        """比軌道エネルギー v^2/2 - GM/r が保存すること。"""
        r = np.linalg.norm(self.position, axis=1)
        v = np.linalg.norm(self.velocity, axis=1)
        energy = 0.5 * v ** 2 - self.gm / r

        drift = np.abs((energy - energy[0]) / energy[0]).max()
        self.assertLess(drift, 1.0e-10, '比エネルギーが保存していない（相対変化 %.3e）' % drift)

    def test_angular_momentum_is_conserved(self):
        """比角運動量 r x v が大きさ・向きとも保存すること。"""
        h = np.cross(self.position, self.velocity)
        h_mag = np.linalg.norm(h, axis=1)

        drift = np.abs((h_mag - h_mag[0]) / h_mag[0]).max()
        self.assertLess(drift, 1.0e-10, '角運動量の大きさが保存していない（相対変化 %.3e）' % drift)

        # 向き（軌道面）も保たれること
        direction = h / h_mag[:, None]
        tilt = np.abs(direction - direction[0]).max()
        self.assertLess(tilt, 1.0e-10, '軌道面が傾いている（%.3e）' % tilt)

    def test_circular_orbit_keeps_its_radius(self):
        """円軌道の半径が一定に保たれること。"""
        r = np.linalg.norm(self.position, axis=1)
        variation = (r.max() - r.min()) / r[0]
        self.assertLess(variation, 1.0e-8, '円軌道の半径が変動している（%.3e）' % variation)

    def test_orbit_returns_to_the_start(self):
        """1 周期後に出発点へ戻ること（周期は解析解 2*pi*sqrt(a^3/GM)）。"""
        r0 = np.linalg.norm(self.position[0])
        period = 2.0 * np.pi * np.sqrt(r0 ** 3 / self.gm)

        # 2 s 刻みなので、周期に最も近いステップを取る
        index = int(round(period / 2.0))
        self.assertLess(index, len(self.position), '積分区間が 1 周期に足りない')

        error = np.linalg.norm(self.position[index] - self.position[0])
        # 1 ステップ分のずれ（周期が刻み幅の整数倍でない）を許容する
        self.assertLess(error, 2.0 * np.linalg.norm(self.velocity[0]) * 2.0,
                        '1 周期後に出発点へ戻っていない（%.3e m）' % error)


class TestRungeKuttaOrder(unittest.TestCase):
    """RK4 が 4 次精度であること。"""

    def test_error_scales_as_fourth_power_of_the_timestep(self):
        gm = gravitational_parameter(circular_orbit_config())

        errors = []
        timesteps = [8.0, 4.0, 2.0]
        for dt in timesteps:
            config = circular_orbit_config()
            _, position, _ = integrate(config, 2000.0, dt)
            # 半径の変動量を誤差の代用にする（厳密解では一定）
            r = np.linalg.norm(position, axis=1)
            errors.append(abs(r[-1] - r[0]))

        # dt を半分にすると誤差は 1/16 になるはず。丸め誤差に埋もれる場合もあるので
        # 「少なくとも 8 倍以上改善する」という緩い条件で確認する
        for coarse, fine in zip(errors[:-1], errors[1:]):
            if coarse < 1.0e-7:
                continue  # すでに丸め誤差レベル
            self.assertGreater(coarse / max(fine, 1.0e-12), 8.0,
                               'RK4 の次数が出ていない: %.3e -> %.3e' % (coarse, fine))


if __name__ == '__main__':
    unittest.main()
