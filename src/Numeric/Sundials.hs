module Numeric.Sundials
  ( -- * The solving function
    solve
    -- * Types
  , OdeProblem(..)
  , Method
  , Tolerances(..)
  , OdeRhsCType
  , SunVector
  , OdeRhs(..)
  , UserData
  , Jacobian
  , ODEOpts(..)
  , SundialsDiagnostics(..)
  , ErrorDiagnostics(..)
  , emptyDiagnostics
  , SundialsSolution(..)
  , EventInfo(..)
  , CrossingDirection(..)
  , EventSpec(..)
  , ARKMethod(..)
  , CVMethod(..)
  ) where

import Numeric.Sundials.Types
import Numeric.Sundials.Common
import Numeric.Sundials.CVode as CV
import Numeric.Sundials.ARKode as ARK
import Control.Monad.IO.Class
import Katip

-- | Solve an ODE system using either ARKode or CVode (depending on what
-- @method@ is instantiated with).
solve
  :: forall m method . (Katip m, Method method)
  => ODEOpts method -- ^ solver options
  -> OdeProblem -- ^ the ODE system to solve
  -> m (Either ErrorDiagnostics SundialsSolution)
solve =
  case methodSolver @method of
    CVode -> solveCommon CV.solveC
    ARKode -> solveCommon ARK.solveC
