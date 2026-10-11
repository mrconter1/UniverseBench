"""Registry of tests, in roadmap order. Each module's docstring describes
its data source and the model method it requires."""
from .clusters import ClusterMasses
from .dwarfs import DwarfSpheroidals
from .milkyway import MilkyWayTest
from .wide_binaries import WideBinaries
from .rotation import RotationExtrapolate, RotationPredict


def all_tests():
    return [RotationPredict(), RotationExtrapolate(), ClusterMasses(), MilkyWayTest(),
            DwarfSpheroidals(), WideBinaries()]
