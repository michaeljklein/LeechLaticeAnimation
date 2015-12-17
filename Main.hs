--import Data.Ratio
import qualified Data.Array.Repa as R
import Control.Monad (liftM)

data Power = Power Double Double
  deriving (Show, Eq, Ord)

dx :: Power -> Power
dx (Power c p) = Power (c*p) (p-1)

diff a b = if a < b then b - a else a - b

matDiff a b = R.sumAllP $ a (R.-^) b

okDiff a b = liftM (((1/100)>).abs) $ matDiff a b

centerAt = 64

eq2 (Z :. y :. x) = if x == y then 1 else 0
eq2c c (Z :. y :. x) = if x == y then c else 0

mId n      = R.fromFunction (ix2 24 24) eq2
mConst n c = R.fromFunction (ix2 24 24) (eq2c c)

-- | gives m^(a/b) within tolerance defined in 'okDiff'. mr = result matrix, mt = this matrix, mp = m0^(some power), nxn input array
--   This is an implementation of a matrix fourier series, centered at 'centerAt'.
matPower m0 a b n = matPower' (Power (centerAt**(a/b)) (a/b)) (mConst n 0) (mId n)
  where
    matPower' (Power c p) mr mp | okDiff mr mr' = mr'
                                | otherwise     = matPower (dx (Power c p)) mr' mp'
                                  where
                                    mp' = m0.mp
                                    mr' = mr + mp'*c
