"""Tests para Hessianos."""

import unittest
import numpy as np
from extremal_packings.configurations import Configuration
from extremal_packings.constraints import build_contact_matrix, rolling_space_basis
from extremal_packings.hessian import (
    build_unconstrained_hessian,
    project_to_roll,
    intrinsic_spectrum
)


class TestBuildUnconstrainedHessian(unittest.TestCase):
    """Tests para build_unconstrained_hessian."""

    def test_dimensions(self):
        """Test dimensiones del Hessiano."""
        coords = np.array([[0, 0], [2, 0], [4, 0]])
        edges = [(0, 1), (1, 2)]
        config = Configuration(coords=coords, edges=edges)
        
        K = build_unconstrained_hessian(config)
        
        # 3 discos -> (6, 6)
        self.assertEqual(K.shape, (6, 6))

    def test_symmetry(self):
        """Test que K sea simétrico."""
        coords = np.array([[0, 0], [2, 0], [1, np.sqrt(3)]])
        edges = [(0, 1), (1, 2), (2, 0)]
        config = Configuration(coords=coords, edges=edges)
        
        K = build_unconstrained_hessian(config)
        
        np.testing.assert_array_almost_equal(K, K.T)

    def test_single_disk(self):
        """Test con un solo disco."""
        coords = np.array([[0, 0]])
        edges = []
        config = Configuration(coords=coords, edges=edges)
        
        K = build_unconstrained_hessian(config)
        
        # Sin contactos -> K debería ser cero
        np.testing.assert_array_almost_equal(K, np.zeros((2, 2)))


class TestProjectToRoll(unittest.TestCase):
    """Tests para project_to_roll."""

    def test_dimensions(self):
        """Test dimensiones de H."""
        coords = np.array([[0, 0], [2, 0], [4, 0]])
        edges = [(0, 1), (1, 2)]
        config = Configuration(coords=coords, edges=edges)
        
        A = build_contact_matrix(config)
        R = rolling_space_basis(A)
        K = build_unconstrained_hessian(config)
        
        H = project_to_roll(K, R)
        
        # H debe ser (d, d) donde d = R.shape[1]
        d = R.shape[1]
        self.assertEqual(H.shape, (d, d))

    def test_symmetry(self):
        """Test que H sea simétrico."""
        coords = np.array([[0, 0], [2, 0], [4, 0]])
        edges = [(0, 1), (1, 2)]
        config = Configuration(coords=coords, edges=edges)
        
        A = build_contact_matrix(config)
        R = rolling_space_basis(A)
        K = build_unconstrained_hessian(config)
        H = project_to_roll(K, R)
        
        np.testing.assert_array_almost_equal(H, H.T)


class TestIntrinsicSpectrum(unittest.TestCase):
    """Tests para intrinsic_spectrum."""

    def test_eigenvalues_ordered(self):
        """Test que autovalores estén ordenados."""
        # Matriz simétrica arbitraria
        H = np.array([
            [3, 1],
            [1, 2]
        ])
        
        eigenvalues = intrinsic_spectrum(H)
        
        # Deben estar ordenados
        self.assertTrue(np.all(eigenvalues[:-1] <= eigenvalues[1:]))

    def test_real_eigenvalues(self):
        """Test que autovalores sean reales."""
        H = np.array([
            [4, 1, 0],
            [1, 3, 1],
            [0, 1, 2]
        ])
        
        eigenvalues = intrinsic_spectrum(H)
        
        # Todos deben ser reales
        self.assertTrue(np.all(np.isreal(eigenvalues)))

    def test_positive_semidefinite(self):
        """Test matriz semidefinida positiva."""
        # Matriz Gramiana -> semidefinida positiva
        A = np.random.randn(5, 3)
        H = A.T @ A
        
        eigenvalues = intrinsic_spectrum(H)
        
        # Todos >= 0
        self.assertTrue(np.all(eigenvalues >= -1e-10))


if __name__ == '__main__':
    unittest.main()
