"""Registry of tests, in roadmap order. Each module's docstring describes
its data source and the model method it requires."""
from .bullet_cluster import BulletCluster
from .clusters import ClusterMasses
from .galaxy_lensing import GalaxyLensing
from .dm_deficient import DMDeficientGalaxies
from .dwarfs import DwarfSpheroidals
from .ellipticals import Ellipticals
from .milkyway import MilkyWayTest
from .wide_binaries import WideBinaries
from .rotation import RotationExtrapolate, RotationPredict


def all_tests():
    return [RotationPredict(), RotationExtrapolate(), ClusterMasses(), MilkyWayTest(),
            DwarfSpheroidals(), WideBinaries(), GalaxyLensing(),
            Ellipticals(), BulletCluster(), DMDeficientGalaxies()]
