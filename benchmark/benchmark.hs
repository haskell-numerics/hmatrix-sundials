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
  { nts :: !Int -- the number of time steps in the output
  , step :: !Double -- the size of a single time step
  , time :: !t -- the total solving time (if known)
  }
  deriving Generic

instance FromNamedRecord (Measurement Time)
instance ToNamedRecord (Measurement Time)
instance DefaultOrdered (Measurement Time)

measure :: Measurement () -> IO (Measurement Time)
measure m@Measurement{..} = do
  let
    rhs _t y =
      let
        mu = 1000
        y1 = y ! 0
        y2 = y ! 1
      in
        fromList [y2, mu * (1 - y1*y1) * y2 - y1]
    y0 = fromList [2,0]
    !times = fromList . take (nts+1) $ iterate (+ step) 0
  time <- (timed . evaluate)
    (odeSolveV BDF Nothing 1e-2 1e-2 rhs y0 times :: Matrix Double)
  return m{time}

main :: IO ()
main = do
  measurements <- sequence $ do
    nts <- [2]
    step <- [1e5, 2e5 .. 2e9]
    let m0 = Measurement { time = (), .. }
    return . unsafeInterleaveIO $ measure m0
  LBS.writeFile "timings.csv" . encodeDefaultOrderedByName . mconcat $ encodeNamedRecord <$> measurements
