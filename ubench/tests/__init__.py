"""Registry of tests, in roadmap order. Each module's docstring describes
its data source and the model method it requires.

STATUS == "settled" tests make up the scored benchmark; "contested" ones are
open questions: reported, never scored, never used to fit constants."""
from .bullet_cluster import BulletCluster
from .clusters import ClusterMasses
from .cosmo import (BAO, BBN, CMB, ClusterCounts, HubbleTension, JWSTEarlyGalaxies, Lithium,
                    LymanAlpha, S8Tension, Supernovae)
from .cosmo_ext import Curvature, GrowthRate, Inflation, Neutrinos, NeutrinoSqueeze, StandardSirens
from .cosmo_open import CosmicDipole, DarkFlow, EarlyBlackHoles
from .dm_deficient import DMDeficientGalaxies
from .dwarfs import DwarfSpheroidals
from .ellipticals import Ellipticals
from .galaxy_lensing import GalaxyLensing
from .gravity import (BinaryPulsars, BlackHoles, Clocks, GravitationalWaves, LabGravity,
                      SolarSystem)
from .milkyway import MilkyWayTest
from .particles import CosmicAntimatter, DirectDetection, IndirectDetection
from .rotation import RotationExtrapolate, RotationPredict
from .wide_binaries import WideBinaries


def all_tests():
    tests = [
        # Stage 1: galaxies and clusters
        RotationPredict(), RotationExtrapolate(), ClusterMasses(), MilkyWayTest(),
        DwarfSpheroidals(), WideBinaries(), GalaxyLensing(), Ellipticals(), BulletCluster(),
        DMDeficientGalaxies(),
        # Stage 2: precision gravity
        SolarSystem(), BinaryPulsars(), GravitationalWaves(), LabGravity(), BlackHoles(), Clocks(),
        # Stage 3: cosmology
        CMB(), BBN(), Lithium(), BAO(), Supernovae(), HubbleTension(), S8Tension(),
        ClusterCounts(), LymanAlpha(), JWSTEarlyGalaxies(),
        # Stage 3b: extended cosmology
        Inflation(), Neutrinos(), NeutrinoSqueeze(), Curvature(), GrowthRate(), CosmicDipole(),
        DarkFlow(), StandardSirens(), EarlyBlackHoles(),
        # Stage 4: particle guard rails
        DirectDetection(), IndirectDetection(), CosmicAntimatter(),
    ]
    return tests
