{-# LANGUAGE FlexibleContexts #-}
--import Data.Ratio
import qualified Data.Array.Repa as R
import Data.Array.Repa.Index (Z, DIM2, ix2)
import Data.Array.Repa.Algorithms.Matrix (mmultS)
import Control.Monad (liftM)

data Power = Power Double Double
  deriving (Show, Eq, Ord)

dx :: Power -> Power
dx (Power c p) = Power (c*p) (p-1)

diff :: (Num a, Ord a) => a -> a -> a
diff a b = if a < b then b - a else a - b

matDiff :: (R.Source r1 Double, R.Source r2 Double) => R.Array r1 DIM2 Double -> R.Array r2 DIM2 Double -> Double
matDiff a b = R.sumAllS $ (R.-^) a b

okDiff :: (R.Source r1 Double, R.Source r2 Double) => R.Array r1 DIM2 Double -> R.Array r2 DIM2 Double -> Bool
okDiff a b = (((0.001)>).abs) $ matDiff a b

centerAt :: Double
centerAt = 64.0

eq2 :: DIM2 -> Double
eq2 (R.Z R.:. y R.:. x) = if x == y then 1.0 else 0.0

eq2c :: Double -> DIM2 -> Double
eq2c c (R.Z R.:. y R.:. x) = if x == y then c else 0.0

mId :: R.Array R.D DIM2 Double
mId = R.fromFunction (ix2 24 24) eq2

mDiag :: Double -> R.Array R.D DIM2 Double
mDiag c = R.fromFunction (ix2 24 24) (eq2c c)

cTimes :: (R.Source r Double) => Double -> R.Array r DIM2 Double -> R.Array R.D DIM2 Double
cTimes c m = R.map (c*) m

-- | gives m^(a/b) within tolerance defined in 'okDiff'. mr = result matrix, mt = this matrix, mp = m0^(some power), nxn input array
--   This is an implementation of a matrix fourier series, centered at 'centerAt'.
matPower :: R.Array R.U DIM2 Double -> Double -> Double -> R.Array R.D DIM2 Double
matPower m0 a b = matPower' (Power (centerAt**(a/b)) (a/b)) (mDiag 0) (R.computeS mId) m0

matPower' :: Power
  -> R.Array R.D DIM2 Double
  -> R.Array R.U DIM2 Double
  -> R.Array R.U DIM2 Double
  -> R.Array R.D DIM2 Double
matPower' (Power c p) mr mp m0  | okDiff mr mr' = mr'
                                | otherwise     = matPower' (dx (Power c p)) mr' mp' m0
  where
    mp' = mmultS m0 mp
    mr' = (R.+^) mr (cTimes c mp')

--Next, need to implement vector-matrix multiplication, hmmm might be able to store a bunch of vectors in a single matrix and just use matrix multiplication...
--Also, need to provide static arrays for leech basis and rotation
--Also, need to interface with (binary?) format for output
