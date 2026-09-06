#!/usr/bin/env python3
"""
姿勢運動の検証（verification）: 厳密解・解析解との比較。

test_attitude.py / test_moment_term.py / test_solver_attitude.py が
「部品が正しいか」「不変条件を満たすか」を見るのに対し、ここでは
**解析的に答えが分かっている問題を production のソルバーに解かせて突き合わせる**。

  1. 収束次数        剛体の等速回転の厳密解に対する次数（RK4 は 4 次、結合系の Euler は 1 次）
  2. 非対称剛体       トルクフリー Euler 方程式の Jacobi 楕円関数による厳密解
  3. 角運動量ベクトル  慣性系で向きまで保存すること（地球自転あり・なしの両方）
  4. 軸対称体の歳差    章動角一定、歳差率 |H|/I_t
  5. 減衰振動         対数減衰率が減衰比 zeta の解析値と一致すること
  6. 重力傾斜 libration 円軌道のピッチ libration 振動数 n*sqrt(3(Ixx-Izz)/Iyy)
  7. 軸対称性         一般の 3 次元姿勢で、力が機体軸と速度の張る面内、
                     モーメントがその面に垂直であること

いずれも数値そのものを見るので、係数・符号・結合の取り違えはここで落ちる。
"""

import contextlib
import copy
import io
import os
import unittest

import numpy as np
from scipy.special import ellipj

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import attitude.attitude as attitude
import coordinate_system.coordinate_system as coordinate_system
import moment_term.moment_term as moment_term
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital


CONFIG_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'config.yml')


def load_config_6dof():
    config = load_config(CONFIG_6DOF)
    config['post_process']['kml']['flag_output'] = False
    config['post_process']['tecplot']['flag_output'] = False
    return config


def vacuum_config(inertia, angular_velocity, attitude_angle=None, rotation_rate=None):
    """
    トルクも空力も無い状態にした config。
    姿勢の運動方程式だけが残るので、剛体の厳密解と直接比べられる。
    """
    config = load_config_6dof()

    config['atmosphere']['kind_atmosphere_model'] = 'constant'
    config['atmosphere'].update({'density': 0.0, 'temperature': 300.0, 'knudsen': 1.0})
    config['satellite']['kind_aerodynamic_model'] = 'constant'
    config['satellite']['drag_coefficient'] = 0.0

    config['attitude'].update({'inertia_tensor': inertia,
                               'flag_moment_aerodynamic': False,
                               'flag_moment_damping': False,
                               'flag_moment_gravity_gradient': False})

    config['initial_settings']['attitude'] = attitude_angle if attitude_angle is not None else [0.0, 0.0, 0.0]
    config['initial_settings']['angular_velocity'] = list(np.array(angular_velocity)*orbital.rad2deg)

    if rotation_rate is not None:
        config['planet']['rotation_rate'] = rotation_rate

    return config


def run_solver(config):
    """production のソルバーを 1 ケース走らせ、姿勢の履歴を返す。"""
    with quiet():
        atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
        aerodynamic_dict = satellite.initial_settings_satellite(config)

        orb = orbital()
        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            orb.initial_settings(config)
        attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)

        iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)

    return {'iteration': iteration, 'coordinate': coordinate_dict, 'velocity': velocity_dict,
            'quaternion': attitude_dict[orbital.KEY_ATTITUDE_QUATERNION],
            'omega': attitude_dict[orbital.KEY_ATTITUDE_OMEGA]}


def pitch_history(config, result):
    """各ステップのピッチ角（deg）を出力ルーチンと同じ手順で作り直す。"""
    angle = []
    for n in range(0, result['iteration']+1):
        polar = coordinate_system.set_angle_polar(config, result['coordinate']['cartesian'][n])
        yaw, pitch, roll = attitude.get_euler_angle(result['quaternion'][n], polar[2], polar[1])
        angle.append(pitch*orbital.rad2deg)
    return np.array(angle)


def zero_crossing(signal, timestep):
    """符号が変わる時刻を線形補間で拾う。"""
    crossing = []
    for n in range(1, len(signal)):
        if signal[n-1]*signal[n] < 0.0:
            fact = signal[n-1]/(signal[n-1] - signal[n])
            crossing.append((n - 1 + fact)*timestep)
    return np.array(crossing)


def rotation_matrix_passive(axis, angle):
    """軸 axis まわりに角 angle だけ座標系を回す受動回転行列（解析解）。"""
    unit = np.array(axis, dtype=float)
    unit = unit/np.linalg.norm(unit)
    cross = np.array([[0.0, -unit[2], unit[1]],
                      [unit[2], 0.0, -unit[0]],
                      [-unit[1], unit[0], 0.0]])
    return np.eye(3) - np.sin(angle)*cross + (1.0 - np.cos(angle))*np.dot(cross, cross)


class TestConvergenceOrder(unittest.TestCase):
    """
    1. 収束次数。

    慣性テンソルが等方でトルクが無ければ角速度は厳密に一定になり、
    姿勢は「角速度まわりの等速回転」という厳密解を持つ。
    その誤差を時間刻みに対して測れば、姿勢側の積分の次数がそのまま出る。
    """

    INERTIA = {'Ixx': 1.0, 'Iyy': 1.0, 'Izz': 1.0}
    OMEGA = np.array([0.3, -0.2, 0.5])
    TIME_MAX = 4.0

    def error_at(self, timestep, scheme):
        config = vacuum_config(self.INERTIA, self.OMEGA, rotation_rate=0.0)
        config['computational_setup']['time_elapsed_maximum'] = self.TIME_MAX
        config['time_integration']['timestep_constant'] = timestep
        config['time_integration']['kind_time_scheme'] = scheme

        result = run_solver(config)

        matrix_initial = attitude.quaternion_to_matrix(result['quaternion'][0])
        matrix_final = attitude.quaternion_to_matrix(result['quaternion'][-1])

        # 初期姿勢から見た相対回転が厳密解と一致するはず
        matrix_relative = np.dot(matrix_final, matrix_initial.T)
        matrix_exact = rotation_matrix_passive(self.OMEGA, np.linalg.norm(self.OMEGA)*self.TIME_MAX)

        return np.linalg.norm(matrix_relative - matrix_exact)

    def test_angular_velocity_stays_exactly_constant(self):
        config = vacuum_config(self.INERTIA, self.OMEGA, rotation_rate=0.0)
        config['computational_setup']['time_elapsed_maximum'] = self.TIME_MAX
        config['time_integration']['timestep_constant'] = 0.1
        result = run_solver(config)
        np.testing.assert_allclose(result['omega'][-1], self.OMEGA, rtol=1.e-14)

    def test_runge_kutta_is_fourth_order(self):
        error = [self.error_at(timestep, 'runge_kutta') for timestep in (0.2, 0.1, 0.05)]
        order = [np.log2(error[n]/error[n+1]) for n in range(0, len(error)-1)]
        for value in order:
            self.assertGreater(value, 3.7)
            self.assertLess(value, 4.3)

    def test_explicit_euler_is_second_order_for_a_pure_rotation(self):
        #
        # オイラー陽解法でも、この問題に限っては 2 次で収束する。
        # dq = dt/2 * Omega q は q と直交するので、1 ステップの更新は
        # 「正しい軸まわりに角度 2 arctan(|omega| dt/2) だけ回す」ことになり、
        # 正規化したあとに残る誤差は角度のずれ
        #   2 arctan(x) - 2x = -2 x^3/3   (x = |omega| dt / 2)
        # だけ、すなわち 1 ステップあたり O(dt^3) になるためである。
        # 正規化をやめると 1 次に落ちるので、これは正規化が効いていることの検査でもある。
        #
        error = [self.error_at(timestep, 'explicit_euler') for timestep in (0.02, 0.01, 0.005)]
        order = [np.log2(error[n]/error[n+1]) for n in range(0, len(error)-1)]
        for value in order:
            self.assertGreater(value, 1.8)
            self.assertLess(value, 2.2)

    def test_coupled_explicit_euler_is_first_order(self):
        # 空力・重力傾斜まで入れた結合系では、素直に 1 次になる
        order = self.coupled_order('explicit_euler', (0.02, 0.01, 0.005))
        for value in order:
            self.assertGreater(value, 0.8)
            self.assertLess(value, 1.3)

    def coupled_order(self, scheme, timestep_list):
        #
        # 空力・重力傾斜まで入れた結合系では厳密解が無いので、
        # 3 つの刻みの差から次数を求める（Richardson）。
        #   p = log2( |y(dt) - y(dt/2)| / |y(dt/2) - y(dt/4)| )
        # 戻り値はクォータニオンと角速度それぞれの次数。
        #
        base = load_config_6dof()
        base['computational_setup']['time_elapsed_maximum'] = 20.0
        base['atmosphere']['kind_atmosphere_model'] = 'constant'
        base['atmosphere'].update({'density': 1.e-8, 'temperature': 1000.0, 'knudsen': 1.0})
        base['satellite']['kind_aerodynamic_model'] = 'constant'
        base['satellite']['drag_coefficient'] = 1.0
        base['attitude']['static_stability_derivative'] = -0.5
        base['attitude']['damping_coefficient'] = {'Clp': -0.1, 'Cmq': -0.5, 'Cnr': -0.5}
        base['initial_settings']['attitude'] = [90.0, 15.0, 0.0]
        base['initial_settings']['angular_velocity'] = [3.0, 0.0, 1.0]
        base['time_integration']['kind_time_scheme'] = scheme

        state = []
        for timestep in timestep_list:
            config = copy.deepcopy(base)
            config['time_integration']['timestep_constant'] = timestep
            result = run_solver(config)
            state.append((np.array(result['quaternion'][-1]), np.array(result['omega'][-1])))

        order = []
        for index in (0, 1):
            difference = [np.linalg.norm(state[n][index] - state[n+1][index]) for n in (0, 1)]
            order.append(np.log2(difference[0]/difference[1]))
        return order

    def test_coupled_runge_kutta_is_fourth_order(self):
        for value in self.coupled_order('runge_kutta', (0.2, 0.1, 0.05)):
            self.assertGreater(value, 3.7)
            self.assertLess(value, 4.3)


class TestTorqueFreeAsymmetricBody(unittest.TestCase):
    """
    2. 非対称剛体のトルクフリー運動。

    Ixx < Iyy < Izz の剛体は Euler 方程式の厳密解を持ち、角速度は
    Jacobi の楕円関数 cn, sn, dn で書ける（Landau & Lifshitz, Mechanics, 37 節）。
    非線形項 omega x (I omega) を含む唯一の閉じた解なので、
    ジャイロ項の符号・係数はここでしか厳密には検査できない。
    """

    INERTIA = {'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0}
    OMEGA = np.array([0.8, 0.0, 0.5])   # cn(0)=1, sn(0)=0, dn(0)=1 に合う初期値
    TIME_MAX = 20.0
    TIMESTEP = 0.02

    @classmethod
    def setUpClass(cls):
        cls.inertia = np.diag([cls.INERTIA['Ixx'], cls.INERTIA['Iyy'], cls.INERTIA['Izz']])

        energy_twice = float(np.dot(cls.OMEGA, np.dot(cls.inertia, cls.OMEGA)))
        momentum_square = float(np.linalg.norm(np.dot(cls.inertia, cls.OMEGA))**2)

        inertia_1, inertia_2, inertia_3 = cls.INERTIA['Ixx'], cls.INERTIA['Iyy'], cls.INERTIA['Izz']
        cls.amplitude = [
            np.sqrt((energy_twice*inertia_3 - momentum_square)/(inertia_1*(inertia_3 - inertia_1))),
            np.sqrt((energy_twice*inertia_3 - momentum_square)/(inertia_2*(inertia_3 - inertia_2))),
            np.sqrt((momentum_square - energy_twice*inertia_1)/(inertia_3*(inertia_3 - inertia_1)))]
        cls.rate = np.sqrt((inertia_3 - inertia_2)*(momentum_square - energy_twice*inertia_1)
                           / (inertia_1*inertia_2*inertia_3))
        cls.modulus = (inertia_2 - inertia_1)*(energy_twice*inertia_3 - momentum_square) \
                    / ((inertia_3 - inertia_2)*(momentum_square - energy_twice*inertia_1))

        config = vacuum_config(cls.INERTIA, cls.OMEGA, rotation_rate=0.0)
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        config['time_integration']['timestep_constant'] = cls.TIMESTEP
        cls.result = run_solver(config)

    @classmethod
    def omega_exact(cls, time):
        sn, cn, dn, phase = ellipj(cls.rate*time, cls.modulus)
        return np.array([cls.amplitude[0]*cn, cls.amplitude[1]*sn, cls.amplitude[2]*dn])

    def test_the_reference_solution_satisfies_eulers_equation(self):
        # 厳密解そのものが正しいことを、Euler 方程式に入れて数値微分で確かめる。
        # 比較の基準を取り違えていたら、この時点で落ちる
        step = 1.e-6
        for time in np.linspace(0.0, self.TIME_MAX, 21):
            omega = self.omega_exact(time)
            derivative = (self.omega_exact(time + step) - self.omega_exact(time - step))/(2.0*step)
            residual = np.dot(self.inertia, derivative) + np.cross(omega, np.dot(self.inertia, omega))
            self.assertLess(np.linalg.norm(residual), 1.e-7)

    def test_initial_condition_matches_the_elliptic_form(self):
        np.testing.assert_allclose(self.omega_exact(0.0), self.OMEGA, atol=1.e-12)

    def test_angular_velocity_matches_the_elliptic_solution(self):
        worst = 0.0
        for n in range(0, self.result['iteration']+1, 10):
            time = float(n)*self.TIMESTEP
            worst = max(worst, np.linalg.norm(np.array(self.result['omega'][n]) - self.omega_exact(time)))
        # 20 秒（楕円関数の周期のおよそ 3 倍）を積分しての差
        self.assertLess(worst, 1.e-9)

    def test_the_solution_is_not_trivial(self):
        # 3 成分とも実際に振れていること（当たり前の解と比べて合格していないことの確認）
        history = np.array(self.result['omega'])
        for index in range(0, 3):
            self.assertGreater(np.ptp(history[:, index]), 0.1)


class TestAngularMomentumInInertialSpace(unittest.TestCase):
    """
    3. 角運動量ベクトルの保存。

    トルクフリーなら角運動量は慣性系で「向きまで」不変である。
    大きさだけでなく向きを見ると、Euler 方程式とクォータニオン積分の
    結合、および地球自転の扱いまでまとめて検査できる。
    """

    INERTIA = {'Ixx': 1.0, 'Iyy': 2.0, 'Izz': 3.0}
    OMEGA = np.array([0.8, 0.3, 0.5])
    TIME_MAX = 20.0
    TIMESTEP = 0.02

    def momentum_history(self, rotation_rate):
        config = vacuum_config(self.INERTIA, self.OMEGA, rotation_rate=rotation_rate)
        config['computational_setup']['time_elapsed_maximum'] = self.TIME_MAX
        config['time_integration']['timestep_constant'] = self.TIMESTEP
        result = run_solver(config)

        inertia = np.diag([self.INERTIA['Ixx'], self.INERTIA['Iyy'], self.INERTIA['Izz']])

        momentum = []
        for n in range(0, result['iteration']+1):
            matrix_be = attitude.quaternion_to_matrix(result['quaternion'][n])
            # 機体軸 -> ECEF
            momentum_ecef = np.dot(matrix_be.T, np.dot(inertia, result['omega'][n]))
            # ECEF -> 慣性系（z 軸まわりに -omega_e t だけ戻す）
            angle = rotation_rate*float(n)*self.TIMESTEP
            matrix_ei = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                                  [np.sin(angle),  np.cos(angle), 0.0],
                                  [0.0, 0.0, 1.0]])
            momentum.append(np.dot(matrix_ei, momentum_ecef))

        return np.array(momentum)

    def test_vector_is_fixed_without_planetary_rotation(self):
        momentum = self.momentum_history(0.0)
        for vector in momentum:
            np.testing.assert_allclose(vector, momentum[0], atol=1.e-10)

    def test_vector_is_fixed_in_inertial_space_with_planetary_rotation(self):
        # 自転を入れると ECEF から見た角運動量は回るが、慣性系では止まっている
        rotation_rate = 7.292115e-5
        momentum = self.momentum_history(rotation_rate)
        for vector in momentum:
            np.testing.assert_allclose(vector, momentum[0], atol=1.e-8)

    def test_the_ecef_frame_really_does_rotate(self):
        # 上の検査が「何もしていないから通った」のではないことの確認
        rotation_rate = 7.292115e-5
        config = vacuum_config(self.INERTIA, self.OMEGA, rotation_rate=rotation_rate)
        config['computational_setup']['time_elapsed_maximum'] = self.TIME_MAX
        config['time_integration']['timestep_constant'] = self.TIMESTEP
        result = run_solver(config)

        inertia = np.diag([self.INERTIA['Ixx'], self.INERTIA['Iyy'], self.INERTIA['Izz']])
        momentum_ecef = []
        for n in (0, result['iteration']):
            matrix_be = attitude.quaternion_to_matrix(result['quaternion'][n])
            momentum_ecef.append(np.dot(matrix_be.T, np.dot(inertia, result['omega'][n])))

        drift = np.linalg.norm(momentum_ecef[-1] - momentum_ecef[0])
        expected = np.linalg.norm(momentum_ecef[0])*rotation_rate*self.TIME_MAX
        self.assertAlmostEqual(drift/expected, 1.0, delta=0.05)


class TestAxisymmetricPrecession(unittest.TestCase):
    """
    4. 軸対称体のトルクフリー歳差。

    慣性系では対称軸が角運動量ベクトルのまわりを一定の章動角で回り、
    その歳差率は |H|/I_t（I_t は横方向の慣性モーメント）になる。
    """

    INERTIA = {'Ixx': 0.75, 'Iyy': 0.50, 'Izz': 0.50}
    OMEGA = np.array([2.0, 0.15, 0.0])
    TIME_MAX = 20.0
    TIMESTEP = 0.005

    @classmethod
    def setUpClass(cls):
        config = vacuum_config(cls.INERTIA, cls.OMEGA, rotation_rate=0.0)
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        config['time_integration']['timestep_constant'] = cls.TIMESTEP
        cls.result = run_solver(config)

        inertia = np.diag([cls.INERTIA['Ixx'], cls.INERTIA['Iyy'], cls.INERTIA['Izz']])
        cls.momentum = np.dot(inertia, cls.OMEGA)
        cls.momentum_magnitude = float(np.linalg.norm(cls.momentum))

        # 慣性系（自転ゼロなので ECEF と同じ）における角運動量と対称軸
        cls.momentum_inertial = []
        cls.axis_inertial = []
        for n in range(0, cls.result['iteration']+1):
            matrix_be = attitude.quaternion_to_matrix(cls.result['quaternion'][n])
            cls.momentum_inertial.append(np.dot(matrix_be.T, np.dot(inertia, cls.result['omega'][n])))
            cls.axis_inertial.append(matrix_be[0, :])   # 機体 x 軸の ECEF 成分
        cls.momentum_inertial = np.array(cls.momentum_inertial)
        cls.axis_inertial = np.array(cls.axis_inertial)

    def test_nutation_angle_is_constant(self):
        unit_momentum = self.momentum_inertial[0]/np.linalg.norm(self.momentum_inertial[0])
        angle = np.arccos(np.clip(np.dot(self.axis_inertial, unit_momentum), -1.0, 1.0))
        self.assertLess(np.ptp(angle), 1.e-9)

        # 章動角の解析値: tan(theta) = I_t w_t / (I_a w_a)
        angle_analytic = np.arctan2(self.INERTIA['Iyy']*np.linalg.norm(self.OMEGA[1:3]),
                                    self.INERTIA['Ixx']*self.OMEGA[0])
        self.assertAlmostEqual(float(angle[0]), float(angle_analytic), places=9)

    def test_precession_rate_matches_the_analytic_value(self):
        # 角運動量に垂直な面に対称軸を投影し、その回転角の進み方を見る
        unit_momentum = self.momentum_inertial[0]/np.linalg.norm(self.momentum_inertial[0])
        basis_1 = np.cross(unit_momentum, [1.0, 0.0, 0.0])
        basis_1 = basis_1/np.linalg.norm(basis_1)
        basis_2 = np.cross(unit_momentum, basis_1)

        phase = np.unwrap(np.arctan2(np.dot(self.axis_inertial, basis_2),
                                     np.dot(self.axis_inertial, basis_1)))
        rate_computed = (phase[-1] - phase[0])/self.TIME_MAX
        rate_analytic = self.momentum_magnitude/self.INERTIA['Iyy']

        self.assertAlmostEqual(abs(rate_computed), rate_analytic, places=6)


class TestDampedOscillation(unittest.TestCase):
    """
    5. 減衰振動。

    定数モデルでは微小迎角のピッチ運動が線形の減衰振動になる:
      I theta'' + c theta' + k theta = 0,
      k = q S L |Cma|,  c = q S L (L/2V) |Cmq|,  zeta = c / (2 sqrt(k I)) .
    振幅の包絡線は exp(-zeta*omega_n*t) で減る。減衰項の大きさはここで初めて数値で検査される。
    （抗力は 0 にして動圧を一定に保ち、減衰が測れるように Cmq を大きめに取ってある）
    """

    DENSITY = 1.e-4
    STABILITY = -0.5
    DAMPING = -20.0
    TIMESTEP = 0.002
    TIME_MAX = 1.2

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        config['time_integration']['timestep_constant'] = cls.TIMESTEP
        config['atmosphere']['kind_atmosphere_model'] = 'constant'
        config['atmosphere'].update({'density': cls.DENSITY, 'temperature': 1000.0, 'knudsen': 1.0})
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['drag_coefficient'] = 0.0     # 減速させず動圧を一定に保つ
        config['attitude']['static_stability_derivative'] = cls.STABILITY
        config['attitude']['damping_coefficient'] = {'Clp': 0.0, 'Cmq': cls.DAMPING, 'Cnr': cls.DAMPING}
        config['attitude']['flag_moment_gravity_gradient'] = False
        config['initial_settings']['attitude'] = [90.0, 2.0, 0.0]
        config['initial_settings']['angular_velocity'] = [0.0, 0.0, 0.0]

        cls.config = config
        cls.result = run_solver(config)
        cls.pitch = pitch_history(config, cls.result)

        velocity = float(np.linalg.norm(cls.result['velocity']['cartesian'][0]))
        area = config['satellite']['characteristic_area']
        length = config['satellite']['characteristic_length']
        inertia = config['attitude']['inertia_tensor']['Iyy']

        dynamic_pressure = 0.5*cls.DENSITY*velocity**2
        cls.stiffness = dynamic_pressure*area*length*abs(cls.STABILITY)
        cls.damping = dynamic_pressure*area*length*length/(2.0*velocity)*abs(cls.DAMPING)
        cls.frequency_natural = np.sqrt(cls.stiffness/inertia)
        cls.ratio_damping = cls.damping/(2.0*np.sqrt(cls.stiffness*inertia))

    def test_the_case_is_underdamped_and_measurable(self):
        self.assertGreater(self.ratio_damping, 0.005)
        self.assertLess(self.ratio_damping, 0.3)

    def test_damped_frequency_matches_the_analytic_value(self):
        crossing = zero_crossing(self.pitch, self.TIMESTEP)
        self.assertGreaterEqual(len(crossing), 4)
        period = 2.0*(crossing[-1] - crossing[0])/float(len(crossing) - 1)

        frequency_analytic = self.frequency_natural*np.sqrt(1.0 - self.ratio_damping**2)
        self.assertAlmostEqual(2.0*np.pi/period/frequency_analytic, 1.0, delta=0.01)

    def test_logarithmic_decrement_matches_the_damping_ratio(self):
        # 極大の絶対値を拾い、ln(amplitude) の傾きから zeta*omega_n を求める
        peak_time = []
        peak_value = []
        for n in range(1, len(self.pitch)-1):
            if abs(self.pitch[n]) >= abs(self.pitch[n-1]) and abs(self.pitch[n]) > abs(self.pitch[n+1]):
                peak_time.append(float(n)*self.TIMESTEP)
                peak_value.append(abs(self.pitch[n]))

        self.assertGreaterEqual(len(peak_value), 4)
        slope = np.polyfit(np.array(peak_time), np.log(np.array(peak_value)), 1)[0]

        decay_analytic = -self.ratio_damping*self.frequency_natural
        self.assertAlmostEqual(slope/decay_analytic, 1.0, delta=0.05)

    def test_amplitude_actually_decays(self):
        self.assertLess(abs(self.pitch).max()*0.8, 2.0)
        self.assertLess(abs(self.pitch[-len(self.pitch)//4:]).max(),
                        0.9*abs(self.pitch[:len(self.pitch)//4]).max())


class TestGravityGradientLibration(unittest.TestCase):
    """
    6. 重力傾斜による libration。

    円軌道で地球指向を保つ剛体の微小ピッチ運動は
      Iyy theta'' + 3 n^2 (Ixx - Izz) theta = 0
    に従い、固有振動数は n sqrt( 3 (Ixx - Izz) / Iyy ) になる（n は軌道角速度）。
    軌道運動・ローカル水平系・クォータニオン・重力傾斜トルクが
    まとめて効く唯一の解析解なので、系全体の検証になる。
    """

    INERTIA = {'Ixx': 2.0, 'Iyy': 1.5, 'Izz': 0.6}
    ALTITUDE = 400.0
    AMPLITUDE = 2.0
    TIMESTEP = 4.0

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()

        # 2 体問題の円軌道にする（J 項と自転を落とす）
        for key in config['planet']['potential_factor']:
            config['planet']['potential_factor'][key] = 0.0
        config['planet']['rotation_rate'] = 0.0

        config['atmosphere']['kind_atmosphere_model'] = 'constant'
        config['atmosphere'].update({'density': 0.0, 'temperature': 300.0, 'knudsen': 1.0})
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['drag_coefficient'] = 0.0

        config['attitude'].update({'inertia_tensor': cls.INERTIA,
                                   'flag_moment_aerodynamic': False,
                                   'flag_moment_damping': False,
                                   'flag_moment_gravity_gradient': True})

        gravitational_parameter = config['planet']['gravitational_constant']*config['planet']['mass']
        radius = config['planet']['radius'] + cls.ALTITUDE*1000.0
        speed = np.sqrt(gravitational_parameter/radius)
        cls.rate_orbit = speed/radius

        cls.frequency_analytic = cls.rate_orbit*np.sqrt(
            3.0*(cls.INERTIA['Ixx'] - cls.INERTIA['Izz'])/cls.INERTIA['Iyy'])

        config['initial_settings']['coordinate'] = [0.0, 0.0, cls.ALTITUDE]
        config['initial_settings']['velocity'] = [speed, 0.0, 0.0]
        # 地球指向を保つには、軌道角速度で機首下げ方向に回り続ける必要がある
        config['initial_settings']['attitude'] = [90.0, cls.AMPLITUDE, 0.0]
        config['initial_settings']['angular_velocity'] = [0.0, -cls.rate_orbit*orbital.rad2deg, 0.0]

        config['computational_setup']['time_elapsed_maximum'] = 1.7*2.0*np.pi/cls.frequency_analytic
        config['time_integration']['timestep_constant'] = cls.TIMESTEP

        cls.config = config
        cls.result = run_solver(config)
        cls.pitch = pitch_history(config, cls.result)

    def test_the_orbit_stays_circular(self):
        radius = [np.linalg.norm(vector) for vector in self.result['coordinate']['cartesian']]
        self.assertLess((max(radius) - min(radius))/radius[0], 1.e-6)

    def test_libration_frequency_matches_the_analytic_value(self):
        crossing = zero_crossing(self.pitch, self.TIMESTEP)
        self.assertGreaterEqual(len(crossing), 3)
        period = 2.0*(crossing[-1] - crossing[0])/float(len(crossing) - 1)

        self.assertAlmostEqual(2.0*np.pi/period/self.frequency_analytic, 1.0, delta=0.01)

    def test_the_libration_stays_small(self):
        # 地球指向姿勢が保たれていること（回り方を取り違えると 90 度級で振れる）
        self.assertLess(np.max(np.abs(self.pitch)), 1.05*self.AMPLITUDE)
        self.assertGreater(np.max(np.abs(self.pitch)), 0.95*self.AMPLITUDE)

    def test_libration_is_slower_than_the_orbit_but_of_the_same_order(self):
        # n sqrt(3 (Ixx-Izz)/Iyy) が軌道角速度の 1.6 倍程度であることの確認
        self.assertAlmostEqual(self.frequency_analytic/self.rate_orbit, np.sqrt(2.8), places=12)


class TestAxisymmetryOfTheAerodynamics(unittest.TestCase):
    """
    7. 軸対称性。

    迎角だけで引く係数表は「回転体」を表している。したがって一般の 3 次元姿勢でも
      - 空力は機体 x 軸と速度が張る面の中にあり、
      - 空力モーメントはその面に垂直である
    のいずれもが厳密に成り立たなければならない。
    面内運動のテストでは検査されない、空力ロール角の回転行列を検証する。
    """

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()
        config['attitude']['center_of_gravity'] = [0.0, 0.0, 0.0]
        cls.config = config
        with quiet():
            cls.aerodynamic_dict = satellite.initial_settings_satellite(config)
            cls.moment, cls.property_dict = moment_term.moment_initialsettings(config)

    def state_at(self, quaternion, velocity):
        return solver.get_aerodynamic_state(
            self.config, self.property_dict, self.aerodynamic_dict, 'fileread',
            np.array([4.0e6, 3.0e6, 3.5e6]), velocity, quaternion, np.zeros(3),
            self.config['satellite']['mass'], self.config['satellite']['characteristic_area'],
            self.config['satellite']['characteristic_length'],
            1.e-7, 1.0, 1.0, 1.0, 0.0,
            0.0, self.moment)

    def test_force_and_moment_respect_the_axis_of_symmetry(self):
        generator = np.random.RandomState(20260906)
        for case in range(0, 20):
            quaternion = attitude.quaternion_normalize(generator.normal(size=4))
            velocity = generator.normal(size=3)
            velocity = 7000.0*velocity/np.linalg.norm(velocity)

            force, moment_total, omega_relative = self.state_at(quaternion, velocity)

            matrix_be = attitude.quaternion_to_matrix(quaternion)
            axis_body = matrix_be[0, :]                     # 機体 x 軸（ECEF 成分）
            normal = np.cross(axis_body, velocity)
            if np.linalg.norm(normal) < 1.e-8:
                continue
            normal = normal/np.linalg.norm(normal)

            # 力は面内: 面の法線と直交する
            force_ecef = np.array(force)
            self.assertLess(abs(np.dot(force_ecef, normal))/np.linalg.norm(force_ecef), 1.e-12)

            # モーメントは面に垂直: 法線と平行
            moment_ecef = np.dot(matrix_be.T, self.moment[1, :])
            if np.linalg.norm(moment_ecef) < 1.e-14:
                continue
            direction = moment_ecef/np.linalg.norm(moment_ecef)
            self.assertLess(np.linalg.norm(np.cross(direction, normal)), 1.e-12)

    def test_the_axial_case_produces_pure_drag(self):
        velocity = np.array([1000.0, -6000.0, 4000.0])
        axis_x = velocity/np.linalg.norm(velocity)
        axis_y = np.cross([0.0, 0.0, 1.0], axis_x)
        axis_y = axis_y/np.linalg.norm(axis_y)
        quaternion = attitude.matrix_to_quaternion(np.array([axis_x, axis_y, np.cross(axis_x, axis_y)]))

        force, moment_total, omega_relative = self.state_at(quaternion, velocity)

        direction = np.array(force)/np.linalg.norm(force)
        np.testing.assert_allclose(direction, -axis_x, atol=1.e-12)
        np.testing.assert_allclose(self.moment[1, :], np.zeros(3), atol=1.e-12)


if __name__ == '__main__':
    unittest.main()
