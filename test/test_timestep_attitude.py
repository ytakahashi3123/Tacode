#!/usr/bin/env python3
"""
姿勢計算の時間刻みについての検証。

姿勢はサブサイクリングを持たず並進と同じ刻みで解くので、「その刻みで振動が
解けているか」がそのまま結果の信頼性になる。solver.check_timestep_attitude は
dt > T/20 で警告を出すが、その T/20 という値の根拠はどこにも測られていなかった。
ここでは 3 つのことを見る。

1. T/20 という基準の妥当性（TestResolutionCriterion）
   解析周期 T の分かっている定数モデルで、dt = T/5, T/10, T/20, T/40 と振り、
   減衰を切ってあるのに振幅が落ちる量（= 数値減衰）を測る。8 周期積んだあとの
   振幅変化は T/5 で -55%、T/10 で -2.8%、T/20 で -0.07% になる。すなわち
   T/20 は「数値減衰が物理のドリフトに埋もれる最初の刻み」で、警告の閾値として
   妥当である。振幅は「状態量から作る不変量」で測る（下記 oscillation_invariant）。

2. 滑らかな空力なら 4 次で収束すること（TestSmoothAerodynamics）
   大気テーブルの内挿を入れたままでも、空力が解析的に滑らかなら RK4 の 4 次が
   出る。すなわち高度方向の内挿は次数を落とさない。

3. 空力テーブルの内挿は C0 でしかないこと（TestTableInterpolation）
   迎角テーブルは 10 deg 刻みの線形内挿なので、ノードで傾きが跳ぶ。迎角は
   振動して同じノードを何度も横切るため、実際のチュートリアル設定では収束次数が
   4 次より落ちる（実測 2〜3 次）。**これは実装の誤りではなくテーブルの性質**で、
   ここで機構そのものを押さえておく。

4. 出荷している刻みが収束していること（TestShippedTimestep）
   チュートリアルの dt = 0.05 s は、dt = 0.025 s と比べて 100 s 後の姿勢が
   1e-5 deg しか違わない。1000 s まで通しても 1e-5 deg 級である（測定値は
   CLAUDE.md に記録した）。
"""

import copy
import os
import unittest

import numpy as np

from context import ROOT_DIR, load_config, quiet

import atmosphere.atmosphere as atmosphere
import attitude.attitude as attitude
import satellite.satellite as satellite
import solver.solver as solver
from orbital.orbital import orbital


CONFIG_6DOF = os.path.join(ROOT_DIR, 'tutorial', 'work_reentry_6dof', 'config.yml')


def load_config_6dof():
    config = load_config(CONFIG_6DOF)
    config['post_process']['kml']['flag_output'] = False
    config['post_process']['tecplot']['flag_output'] = False
    return config


def run_solver(config):
    """production のソルバーを 1 ケース走らせ、履歴を返す。"""
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


def final_state(config, timestep):
    config = copy.deepcopy(config)
    config['time_integration']['timestep_constant'] = timestep
    result = run_solver(config)
    return np.array(result['quaternion'][-1]), np.array(result['omega'][-1])


def attitude_difference(quaternion, quaternion_reference):
    """2 つの姿勢の隔たりを回転角（deg）で測る。"""
    matrix_relative = np.dot(attitude.quaternion_to_matrix(quaternion),
                             attitude.quaternion_to_matrix(quaternion_reference).T)
    cosine = 0.5*(np.trace(matrix_relative) - 1.0)
    return np.arccos(min(1.0, max(-1.0, cosine)))*orbital.rad2deg


def convergence_order(config, timestep_list, index):
    """
    3 つの刻みの差から収束次数を求める（Richardson）。
      p = log2( |y(dt) - y(dt/2)| / |y(dt/2) - y(dt/4)| )
    index は 0 でクォータニオン、1 で角速度。
    """
    state = [final_state(config, timestep) for timestep in timestep_list]
    difference = [np.linalg.norm(state[n][index] - state[n+1][index])
                  for n in range(0, len(state)-1)]
    return [np.log2(difference[n]/difference[n+1]) for n in range(0, len(difference)-1)]


class TestResolutionCriterion(unittest.TestCase):
    """
    T/20 という判定基準の妥当性。

    定数モデルで減衰も重力傾斜も切ってあるので、微小迎角のピッチ振動は
      T = 2 pi sqrt( Iyy / (q S L |Cm_alpha|) )
    の等振幅振動になる。振幅が落ちればそれは数値減衰である。
    """

    DENSITY = 1.e-8
    STABILITY = -0.5
    ANGLE_INITIAL = 2.0
    NUMBER_PERIOD = 8
    DIVISOR = (5, 10, 20, 40)

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()
        config['atmosphere']['kind_atmosphere_model'] = 'constant'
        config['atmosphere'].update({'density': cls.DENSITY, 'temperature': 1000.0,
                                     'knudsen': 1.0})
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['drag_coefficient'] = 1.0
        config['attitude']['static_stability_derivative'] = cls.STABILITY
        config['attitude']['flag_moment_damping'] = False
        config['attitude']['flag_moment_gravity_gradient'] = False
        config['initial_settings']['attitude'] = [90.0, cls.ANGLE_INITIAL, 0.0]
        # 地球自転を止めて、慣性系基準と ECEF 基準の差を持ち込まない
        config['planet']['rotation_rate'] = 0.0

        dynamic_pressure = 0.50*cls.DENSITY*config['initial_settings']['velocity'][0]**2
        stiffness = abs(dynamic_pressure*config['satellite']['characteristic_area']
                        * config['satellite']['characteristic_length']*cls.STABILITY)
        cls.omega_natural = np.sqrt(stiffness/config['attitude']['inertia_tensor']['Iyy'])
        cls.period = 2.0*np.pi/cls.omega_natural

        config['computational_setup']['time_elapsed_maximum'] = cls.NUMBER_PERIOD*cls.period
        cls.config = config

        cls.measured = {}
        for divisor in cls.DIVISOR:
            cls.measured[divisor] = cls.measure(cls.period/divisor, divisor)

    @classmethod
    def measure(cls, timestep, divisor):
        """1 つの刻みで走らせ、周期と振幅の変化を測る。"""
        config = copy.deepcopy(cls.config)
        config['time_integration']['timestep_constant'] = timestep
        result = run_solver(config)

        angle = []
        invariant = []
        for n in range(0, result['iteration']+1):
            quaternion = result['quaternion'][n]
            matrix_be = attitude.quaternion_to_matrix(quaternion)
            velocity_body = np.dot(matrix_be, result['velocity']['cartesian'][n])
            alpha, beta, alpha_total, phi_aero = attitude.get_aerodynamic_angle(velocity_body)
            omega_relative = attitude.get_omega_relative(result['omega'][n], 0.0, matrix_be)
            angle.append(alpha)
            #
            # 振幅の指標は「状態量から作る不変量」にする
            #   A = sqrt( alpha^2 + (q/omega_n)^2 )
            # 等振幅の調和振動なら常に振幅そのものに等しい。ピークの値を拾うと
            # 粗い刻みでは頂点を跨いで標本化してしまい、数値減衰と区別できない。
            #
            invariant.append(np.sqrt(alpha**2 + (omega_relative[1]/cls.omega_natural)**2))

        angle = np.array(angle)
        invariant = np.array(invariant)*orbital.rad2deg

        # 周期は符号反転の間隔から。半周期ごとに 1 回横切る
        crossing = [(n - 1 + angle[n-1]/(angle[n-1] - angle[n]))*timestep
                    for n in range(1, len(angle)) if angle[n-1]*angle[n] < 0.0]
        period = 2.0*np.mean(np.diff(crossing)) if len(crossing) >= 3 else float('nan')

        # 最初と最後の 1 周期で不変量を平均して比べる
        amplitude_first = float(np.mean(invariant[0:divisor+1]))
        amplitude_last = float(np.mean(invariant[-divisor-1:]))

        return {'period_error': (period - cls.period)/cls.period,
                'amplitude_change': amplitude_last/amplitude_first - 1.0}

    def test_the_oscillation_is_actually_resolved_at_all(self):
        # 8 周期を積んで周期が測れていること（テストが空回りしていないこと）
        for divisor in self.DIVISOR:
            self.assertFalse(np.isnan(self.measured[divisor]['period_error']),
                             'dt = T/{} で周期が測れていない'.format(divisor))

    def test_a_coarse_timestep_damps_the_oscillation_spuriously(self):
        # 減衰を切ってあるのに、T/5 では 8 周期で振幅が 3 割以上落ちる
        self.assertLess(self.measured[5]['amplitude_change'], -0.30)
        # T/10 でもまだ 0.5% 以上落ちる
        self.assertLess(self.measured[10]['amplitude_change'], -5.e-3)

    def test_the_criterion_of_one_twentieth_removes_the_numerical_damping(self):
        # T/20 まで細かくすると、8 周期での振幅変化は 0.2% を下回る
        self.assertLess(abs(self.measured[20]['amplitude_change']), 2.e-3)
        # さらに半分にしても振幅の挙動は変わらない（収束している）
        self.assertLess(abs(self.measured[40]['amplitude_change']), 2.e-3)

    def test_the_numerical_damping_shrinks_as_the_timestep_shrinks(self):
        change = [abs(self.measured[divisor]['amplitude_change']) for divisor in (5, 10, 20)]
        self.assertLess(change[1], change[0])
        self.assertLess(change[2], change[1])

    def test_the_period_at_one_twentieth_agrees_with_the_converged_value(self):
        # T/20 の周期は T/40 の周期と 0.1% 以内で一致する。
        # なお両者に共通して残る -3e-4 のずれは離散化ではなく物理側
        # （抗力で減速して動圧が下がる・有限振幅の効果）である
        self.assertLess(abs(self.measured[20]['period_error']
                            - self.measured[40]['period_error']), 1.e-3)
        # 解析周期そのものからのずれも 0.1% 以内
        self.assertLess(abs(self.measured[20]['period_error']), 1.e-3)

    def test_the_solver_warns_about_a_coarse_timestep(self):
        # check_timestep_attitude が T/20 より粗い刻みで実際に警告すること
        config = copy.deepcopy(self.config)
        config['computational_setup']['time_elapsed_maximum'] = 2.0*self.period
        config['time_integration']['timestep_constant'] = self.period/5.0
        with quiet() as buffer:
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)
        self.assertIn('coarse for the attitude oscillation', buffer.getvalue())


class TestSmoothAerodynamics(unittest.TestCase):
    """
    空力が解析的に滑らかなら、大気テーブルの内挿を入れたままでも 4 次で収束する。
    すなわち高度方向の内挿は収束次数を落とさない（高度は単調に進むので、
    テーブルの節点を横切る回数が少ない）。
    """

    TIME_MAX = 100.0
    TIMESTEP = (0.4, 0.2, 0.1)

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        # 大気は fileread のまま。空力だけ定数モデル（滑らか）にする
        config['satellite']['kind_aerodynamic_model'] = 'constant'
        config['satellite']['drag_coefficient'] = 1.0
        config['attitude']['static_stability_derivative'] = -0.5
        cls.config = config

    def test_runge_kutta_is_fourth_order_for_the_quaternion(self):
        for order in convergence_order(self.config, self.TIMESTEP, 0):
            self.assertGreater(order, 3.7)
            self.assertLess(order, 4.3)

    def test_runge_kutta_is_fourth_order_for_the_angular_velocity(self):
        for order in convergence_order(self.config, self.TIMESTEP, 1):
            self.assertGreater(order, 3.7)
            self.assertLess(order, 4.3)


class TestTableInterpolation(unittest.TestCase):
    """
    空力テーブルの内挿は C0 でしかない。

    迎角の節点は 10 deg 刻みで、そこで係数の傾きが跳ぶ。迎角は振動して同じ節点を
    何度も横切るので、テーブルを使った計算の収束次数は 4 次より落ちる（実測 2〜3 次）。
    実装の誤りではなくテーブルの性質なので、次数そのものではなく機構を押さえる。
    次数を回復させたければ、テーブルを細かくするか内挿を高次にする必要がある。
    """

    KNUDSEN = 82.7
    NODE = (10.0, 20.0)
    DELTA = 1.e-3

    @classmethod
    def setUpClass(cls):
        with quiet():
            cls.aerodynamic_dict = satellite.initial_settings_satellite(load_config_6dof())

    def moment_slope(self, angle_lower, angle_upper):
        """ピッチングモーメント係数の差分商。"""
        coefficient = []
        for angle in (angle_lower, angle_upper):
            force, moment = satellite.get_aerodynamic_coefficient_attitude(
                self.KNUDSEN, angle, self.aerodynamic_dict)
            coefficient.append(moment[1])
        return (coefficient[1] - coefficient[0])/(angle_upper - angle_lower)

    def test_the_coefficient_is_continuous_at_a_node(self):
        for node in self.NODE:
            force_lower, moment_lower = satellite.get_aerodynamic_coefficient_attitude(
                self.KNUDSEN, node - self.DELTA, self.aerodynamic_dict)
            force_upper, moment_upper = satellite.get_aerodynamic_coefficient_attitude(
                self.KNUDSEN, node + self.DELTA, self.aerodynamic_dict)
            self.assertAlmostEqual(moment_lower[1], moment_upper[1], places=4)

    def test_the_slope_jumps_at_a_node(self):
        # 節点の左右で片側傾きが 10% 以上違う。これが 4 次収束を崩す原因
        for node in self.NODE:
            slope_lower = self.moment_slope(node - self.DELTA, node)
            slope_upper = self.moment_slope(node, node + self.DELTA)
            self.assertGreater(abs(slope_upper - slope_lower)/abs(slope_lower), 0.10,
                               'AOA = {} deg で傾きが跳んでいない'.format(node))

    def test_the_slope_is_smooth_away_from_a_node(self):
        # 節点から離れた区間の中では線形なので、傾きは一致する
        slope_lower = self.moment_slope(12.0, 13.0)
        slope_upper = self.moment_slope(16.0, 17.0)
        self.assertAlmostEqual(slope_lower, slope_upper, places=10)


class TestShippedTimestep(unittest.TestCase):
    """
    チュートリアルが出荷している dt = 0.05 s で、姿勢が収束していること。

    ここだけはテーブルも含めたチュートリアルそのままの設定で見る。窓は 100 s
    （全長 1000 s の測定値は CLAUDE.md に記録した）。
    """

    TIME_MAX = 100.0
    TIMESTEP_SHIPPED = 0.05
    TIMESTEP_REFERENCE = 0.025
    TOLERANCE_ANGLE = 1.e-4      # deg
    TOLERANCE_OMEGA = 1.e-6      # rad/s

    @classmethod
    def setUpClass(cls):
        config = load_config_6dof()
        config['computational_setup']['time_elapsed_maximum'] = cls.TIME_MAX
        cls.config = config
        cls.reference = final_state(config, cls.TIMESTEP_REFERENCE)
        cls.shipped = final_state(config, cls.TIMESTEP_SHIPPED)
        cls.coarse = final_state(config, 2.0*cls.TIMESTEP_SHIPPED)

    def test_the_shipped_timestep_is_converged(self):
        angle = attitude_difference(self.shipped[0], self.reference[0])
        self.assertLess(angle, self.TOLERANCE_ANGLE,
                        '出荷している dt でも姿勢が {:.3e} deg ずれている'.format(angle))
        self.assertLess(float(np.linalg.norm(self.shipped[1] - self.reference[1])),
                        self.TOLERANCE_OMEGA)

    def test_halving_the_timestep_improves_the_answer(self):
        # dt を半分にした側のほうが基準に近いこと（収束していることの確認）
        self.assertLess(attitude_difference(self.shipped[0], self.reference[0]),
                        attitude_difference(self.coarse[0], self.reference[0]))

    def test_the_solver_does_not_warn_at_the_shipped_timestep(self):
        # dt = 0.05 s はこの窓のあいだ T/20 を割らないので警告は出ない
        config = copy.deepcopy(self.config)
        config['computational_setup']['time_elapsed_maximum'] = 20.0
        config['time_integration']['timestep_constant'] = self.TIMESTEP_SHIPPED
        with quiet() as buffer:
            atmosphere_dict = atmosphere.initial_settings_atmosphere(config)
            aerodynamic_dict = satellite.initial_settings_satellite(config)
            orb = orbital()
            iteration, time_elapsed, coordinate_dict, velocity_dict, trajectory_dict = \
                orb.initial_settings(config)
            attitude_dict = orb.initial_settings_attitude(config, coordinate_dict, velocity_dict)
            solver.solve_equation_motion(config, iteration, time_elapsed,
                                         coordinate_dict, velocity_dict, trajectory_dict,
                                         atmosphere_dict, aerodynamic_dict, attitude_dict)
        self.assertNotIn('coarse for the attitude oscillation', buffer.getvalue())


if __name__ == '__main__':
    unittest.main()
