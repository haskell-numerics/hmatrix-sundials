module Numeric.Sundials
  ( -- * Solving functions
    solveCV
  , solveARK
    -- * Types
  , OdeProblem(..)
  , Tolerances(..)
  , OdeRhsCType
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

solveCV
  :: Katip m
  => ODEOpts CVMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCV = solveCommon CV.solveC

solveARK
  :: Katip m
  => ODEOpts ARKMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveARK = solveCommon ARK.solveC
