import Prelude hiding ((<>))
import Numeric.LinearAlgebra hiding (size, step)
import qualified System.Clock as Clock
import Data.Csv (FromNamedRecord, ToNamedRecord, DefaultOrdered)
import Data.Csv.Incremental
import GHC.Generics (Generic)
import qualified Data.ByteString.Lazy as LBS
import Numeric.Sundials.CVode.ODE
import Control.Exception
import System.IO.Unsafe

type Time = Double

getTime :: IO Time
getTime = do
  t <- Clock.getTime Clock.Monotonic
  let ns = realToFrac $ Clock.toNanoSecs t
  return $ ns / 10 ^ (9 :: Int)

timed :: IO a -> IO Time
timed t = do
  start <- getTime
  !_    <- t
  end   <- getTime
  return (end-start)

data Measurement t = Measurement
  { size :: !Int -- the number of variables/equations in the ODE system
  , nts :: !Int -- the number of time steps in the output
  , step :: !Double -- the size of a single time step
  , time :: !t -- the total solving time (if known)
  }
  deriving Generic

instance FromNamedRecord (Measurement Time)
instance ToNamedRecord (Measurement Time)
instance DefaultOrdered (Measurement Time)

-- | Generate a random Gaussian vector
randVec :: Int -> IO (Vector Double)
randVec n = (! 0) <$> (randn 1 n)

-- | Generate a random negative definite matrix
randNeg :: Int -> IO (Matrix Double)
randNeg n = do
  q <- randn n n
  v <- abs <$> randVec n
  return (- (tr q <> diag v <> q))

measure :: Matrix Double -> Vector Double -> Measurement () -> IO (Measurement Time)
measure odeMatrix y0 m@Measurement{..} = do
  let
    rhs _t y = odeMatrix #> y
    !times = fromList . take (nts+1) $ iterate (+ step) 0
  time <- (timed . evaluate)
    (odeSolveV BDF Nothing 1e-2 1e-2 rhs y0 times :: Matrix Double)
  return m{time}

main :: IO ()
main = do
  measurements <- sequence $ do
    size <- [1 .. 40]
    let
      -- Of course this should be done using some sort of streaming, but
      -- cassava's incremental encoding seems to rely on lazy io.
      !y0 = unsafePerformIO $ randVec size
      !mx = unsafePerformIO $ randNeg size
    nts <- [100, 200 .. 1000]
    step <- [1e3, 2e3 .. 5e6]
    let m0 = Measurement { time = (), .. }
    return . unsafeInterleaveIO $ measure mx y0 m0
  LBS.writeFile "timings.csv" . encodeDefaultOrderedByName . mconcat $ encodeNamedRecord <$> measurements
