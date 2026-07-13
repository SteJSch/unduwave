"""
Unduwave init statements
"""

from unduwave.unduwave_incl import *
from unduwave.constants import *
import unduwave.constants as constants
from unduwave.api.wave_api_root import wave_api as wave
from unduwave.api.undu_api_root import undu_api as undu
from .attribute_classes.attributes import *
import unduwave.quantities.quantities as quantities
from .analytical_module.ana_undu import bfield as bfield

import unduwave.analytical_module.ana_undu.grids as grids
import unduwave.analytical_module.ana_undu.kicks as kicks

from .analytical_module.ana_undu import undulator as undulator
from .wave_modules.wave_parameters import *
from .wave_modules.wave_postprocess import wave_postprocess
from .wave_modules.wave_prepare import wave_prepare
from .wave_modules.wave_control import wave_control
from .wave_modules.wave_results import wave_results
# import unduwave.helpers.bfield_helpers as bfield
# import unduwave.helpers.bfield_helpers as bfield

import unduwave.helpers.file_folder_helpers as ff_h
import unduwave.helpers.numerical_helpers as n_h
# from unduwave.helpers.file_folder_helpers import loadBessyIIundulatorList
import unduwave.undu_modules.undu_parameters as undu_parameters
from .undu_modules.undu_postprocess import undu_postprocess
from .undu_modules.undu_prepare import undu_prepare
from .undu_modules.undu_control import undu_control
from .undu_modules.undu_results import undu_results
# from .undu_list.undu_list_loading import undu_list_loading
import unduwave.undu_list.undu_list_loading as undu_list_loading

import unduwave.undu_modules.undu_blocks as undu_blocks
import unduwave.undu_modules.undu_coils as undu_coils

import unduwave.undu_modules.undu_magneticObjectGeometries as magneticObjectGeometries
import unduwave.undu_modules.undu_undulatorComponents as undulatorComponents

# import unduwave.undu_modules.undu_representation as undu_representation

from .analytical_module.ana_wave import ana_eb_fields as ana_eb
from .analytical_module.ana_undu import undulator as undulator

# predefined undulators

# import unduwave.undulatorsBessyII.cpmu20 as cpmu20
