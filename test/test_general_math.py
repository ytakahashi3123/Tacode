#!/usr/bin/env python3
"""
`general.py` に置いた小さなベクトル演算のテスト。

`cross_product` と `vector_norm` は numpy の `np.cross` / `np.linalg.norm` を
**速度のためだけに**置き換えたものなので、検査することは 1 つしかない:
**値がビット単位で一致すること**。

なぜ置き換えたか（2026-09-16、profile に基づく）:

  - `np.cross` は任意次元・任意軸に対応するので `moveaxis` と
    `normalize_axis_tuple` を通る。3 成分では呼び出しの手間が計算そのものを覆い、
    6 自由度計算の 20 % がここだった（実測 7.3 us -> 0.7 us）
  - `np.linalg.norm` は 1 次元・ord 既定のとき結局 `sqrt(dot(x,x))` を計算している
    （実測 0.82 us -> 0.62 us）

どちらも演算の順序を変えていないので一致する。ここが崩れると、チュートリアルと
validation の出力がバイト一致しなくなる（それが高速化の合格条件だった）。
"""

import unittest

import numpy as np

from context import quiet  # noqa: F401  (import path を通すために読む)

from general.general import cross_product, vector_norm


class TestCrossProduct(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(20260916)

    def test_it_matches_numpy_bit_for_bit(self):
        for _ in range(0, 2000):
            u = self.rng.normal(size=3)*10.0**self.rng.integers(-12, 12)
            v = self.rng.normal(size=3)*10.0**self.rng.integers(-12, 12)
            np.testing.assert_array_equal(cross_product(u, v), np.cross(u, v))

    def test_the_basis_vectors_come_out_right_handed(self):
        ex = np.array([1.0, 0.0, 0.0])
        ey = np.array([0.0, 1.0, 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        np.testing.assert_array_equal(cross_product(ex, ey), ez)
        np.testing.assert_array_equal(cross_product(ey, ez), ex)
        np.testing.assert_array_equal(cross_product(ez, ex), ey)

    def test_it_takes_a_list_as_numpy_does(self):
        # 機体の重心オフセットのように、config 由来のリストが渡ることがある
        np.testing.assert_array_equal(cross_product([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
                                      np.cross([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]))


class TestVectorNorm(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(20260916)

    def test_it_matches_numpy_bit_for_bit(self):
        for length in (3, 4):
            # 3 成分は位置・速度、4 成分はクォータニオン
            for _ in range(0, 2000):
                vector = self.rng.normal(size=length)*10.0**self.rng.integers(-12, 12)
                self.assertEqual(vector_norm(vector), np.linalg.norm(vector))

    def test_a_zero_vector_is_zero(self):
        self.assertEqual(vector_norm(np.zeros(3)), 0.0)


if __name__ == '__main__':
    unittest.main()
