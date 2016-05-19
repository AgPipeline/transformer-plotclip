#!/usr/bin/env python

'''
EnvironmentalLoggerAnalyser.py

----------------------------------------------------------------------------------------
This module will read data generated by Environmental Sensor and convert to netCDF file
----------------------------------------------------------------------------------------
Prerequisite:
1. Python (2.7+ recommended)
2. netCDF4 module for Python (and its dependencies)
3. numpy (For array calculations, make sure the numpy has the same Python verison as other modules)
----------------------------------------------------------------------------------------

Usage: Give full path to EnvironmentalLoggerAnalyser.py, or place it in PYTHONPATH, then:

python EnvironmentalLoggerAnalyser.py drc_in drc_out # Process all files in drc_in
python EnvironmentalLoggerAnalyser.py  fl_in drc_out # Process only fl_in
where drc_in is input directory, drc_out is output directory, fl_in is input file
Input  filenames must have '.json' extension
Output filenames are replace '.json' with '.nc'

UCI test:
python ${HOME}/terraref/computing-pipeline/scripts/hyperspectral/EnvironmentalLoggerAnalyser.py ${DATA}/terraref/environmentlogger_test.json ${DATA}/terraref

UCI production:
python ${HOME}/terraref/computing-pipeline/scripts/hyperspectral/EnvironmentalLoggerAnalyser.py ${DATA}/terraref/EnvironmentLogger/2016-04-07/2016-04-07_12-00-07_enviromentlogger.json ~/rgr

Roger production:
module add gdal-stack-2.7.10
python ${HOME}/terraref/computing-pipeline/scripts/hyperspectral/EnvironmentalLoggerAnalyser.py /projects/arpae/terraref/raw_data/ua-mac/EnvironmentLogger/2016-04-07/2016-04-07_12-00-07_enviromentlogger.json ~/rgr

EnvironmentalLoggerAnalyser.py takes the first argument as the input folder (containing JSON files,
but it can also be one single file) and the second argument as the output folder (output netCDF files go here).
If the output folder does not exist, EnvironmentalLoggerAnalyser.py creates it.

----------------------------------------------------------------------------------------
20160429: Output JSON file is now completely composed by variables
          2D spectrometer variables (wavelength and spectrum) are available in the exported file
20160503: Add chunksizes parameters for time, which significantly reduces the processing time (and the file size)
          Add timestamps and commandLine the user used for each exported file
          Remind the user currently the script is dealing with which file
20160508: 1. "timestamp" array are double-precison, with unit of days offset from UNIX base time (Gregorian Calender)
          2. Remove unnecessary I/O (in wavelength, because we only need the wavelengths in the first set of readings)
20160509: 1. Retrieve the adjustment Professor Zender made in version a38ca7d, May 4th
          2. Rename wavelength variable and dimension as "wvl_lgr" to avoid name collision with Hyperspectral wavelengths
20160511: Add variables and calculations from terraref.nco, including wvl_dlt and flx_sns
          Add new array flx_spc_dwn (downwelling spectral flux) to the netCDF file
20160512: Recalculate downwellingSpectralFlux, save Flux sensitivity in SI
20160517: Implement variable wavelength
20160518: 1. Add numpy module, now the array calculation will be done by numpy to improve efficiency
          2. Add Downwelling Flux (the previous one is recognized and renamed as Downwelling spectral flux)
20160519: 1. Recalculate and double check the method used for calculating downwelling spectral flux
          2. Reinstate the integration time and area (based on the discussion about the dimension of the flux sensitivity)
          3. Clean up based on Professor Zender's adjustment

TODO:
1. Convert all units to SI
----------------------------------------------------------------------------------------
Note:
If you need a different base time, it is named "_UNIX_BASETIME" and located at the
beginning of the script as a global variable. You could simply change the parameters
as they are named.

----------------------------------------------------------------------------------------
Thanks for the advice from Professor Zender and testing data from Mr. Maloney
----------------------------------------------------------------------------------------
'''
import numpy as np
import json
import time
import sys
import os
from datetime import date, datetime
from netCDF4 import Dataset

_UNIT_DICTIONARY = {u'm': 'meter', u"hPa": "hectopascal", u"DegCelsius": "celsius",
                    u's': 'second', u'm/s': 'meter second-1', u"mm/h": 'millimeter hour-1',
                    u"relHumPerCent": "percent", u"?mol/(m^2*s)": "micromole meter-2 second-1", u"umol/(m^2*s)": "micromole meter-2 second-1",
                    u'kilo Lux': 'kilolux', u'degrees': 'degrees', u'?s': 'microsecond', u'us': 'microsecond', '': ''}
_NAMES = {'sensor par': 'Sensor Photosynthetically Active Radiation'}
_UNIX_BASETIME = date(year=1970, month=1, day=1)
_FLX_SNS = \
    [2.14905162e-3, 2.14905162e-3, 2.13191329e-3, 2.09958721e-3, 2.07255014e-3, 2.04042869e-3, 2.01065201e-3, 1.98247712e-3, 1.95695595e-3, 1.93084270e-3,
     1.90570597e-3, 1.87482001e-3, 1.85556117e-3, 1.83111292e-3, 1.80493827e-3, 1.78204243e-3, 1.76011709e-3, 1.74298970e-3, 1.72493868e-3, 1.70343113e-3,
     1.68553722e-3, 1.66631622e-3, 1.65443041e-3, 1.63642311e-3, 1.61846215e-3, 1.60545573e-3, 1.59131286e-3, 1.57713520e-3, 1.56105571e-3, 1.54524417e-3,
     1.53346351e-3, 1.51857527e-3, 1.50523348e-3, 1.49241475e-3, 1.47961893e-3, 1.46809459e-3, 1.45371019e-3, 1.44115770e-3, 1.42656173e-3, 1.41449893e-3,
     1.40218325e-3, 1.38727829e-3, 1.37516400e-3, 1.36250816e-3, 1.34930545e-3, 1.33658553e-3, 1.32360948e-3, 1.31037104e-3, 1.29485640e-3, 1.28218362e-3,
     1.26643446e-3, 1.25134745e-3, 1.23774104e-3, 1.22154062e-3, 1.20610594e-3, 1.19103569e-3, 1.17634593e-3, 1.16196432e-3, 1.14706609e-3, 1.13382928e-3,
     1.11918389e-3, 1.10556178e-3, 1.09188557e-3, 1.07918790e-3, 1.06742601e-3, 1.05718118e-3, 1.04694473e-3, 1.03659582e-3, 1.02873900e-3, 1.02288609e-3,
     1.01628109e-3, 1.01158157e-3, 1.00762812e-3, 1.00491506e-3, 9.99889414e-4, 9.95638736e-4, 9.91223764e-4, 9.86536043e-4, 9.82378851e-4, 9.75992983e-4,
     9.68055921e-4, 9.60360601e-4, 9.52229616e-4, 9.44288731e-4, 9.35551583e-4, 9.27778978e-4, 9.19234295e-4, 9.10434775e-4, 9.02072493e-4, 8.94082691e-4,
     8.86054070e-4, 8.80176269e-4, 8.73884929e-4, 8.67741565e-4, 8.61294713e-4, 8.55225782e-4, 8.50034434e-4, 8.44592095e-4, 8.39026291e-4, 8.32991645e-4,
     8.27232672e-4, 8.20422797e-4, 8.12627476e-4, 8.03944155e-4, 7.93686106e-4, 7.82697718e-4, 7.71489658e-4, 7.58484857e-4, 7.47988099e-4, 7.37636214e-4,
     7.28366631e-4, 7.19162810e-4, 7.11125553e-4, 7.03209565e-4, 6.94909585e-4, 6.86758690e-4, 6.78942712e-4, 6.70144668e-4, 6.62797996e-4, 6.54343779e-4,
     6.46170305e-4, 6.38872127e-4, 6.30959379e-4, 6.25264821e-4, 6.19260059e-4, 6.14223740e-4, 6.09600378e-4, 6.04109797e-4, 5.99684341e-4, 5.94762462e-4,
     5.89358438e-4, 5.84542583e-4, 5.79374467e-4, 5.76250247e-4, 5.72012513e-4, 5.68631187e-4, 5.65821491e-4, 5.62718648e-4, 5.59116033e-4, 5.55411029e-4,
     5.51801247e-4, 5.47598043e-4, 5.43247327e-4, 5.37762461e-4, 5.33460619e-4, 5.28878572e-4, 5.23689536e-4, 5.19266009e-4, 5.14458405e-4, 5.10751276e-4,
     5.07003427e-4, 5.03165882e-4, 4.98966826e-4, 4.96001001e-4, 4.93715941e-4, 4.90249980e-4, 4.86771807e-4, 4.83624709e-4, 4.80476687e-4, 4.76929878e-4,
     4.74001018e-4, 4.70608674e-4, 4.66770209e-4, 4.64581056e-4, 4.61759438e-4, 4.59189259e-4, 4.57332260e-4, 4.55348303e-4, 4.53789335e-4, 4.51997104e-4,
     4.51229308e-4, 4.49813718e-4, 4.49043888e-4, 4.48243809e-4, 4.47227105e-4, 4.45901662e-4, 4.44320515e-4, 4.42138136e-4, 4.40911947e-4, 4.38857020e-4,
     4.36114912e-4, 4.32791420e-4, 4.29184734e-4, 4.26039764e-4, 4.23760025e-4, 4.20638491e-4, 4.17739943e-4, 4.15334137e-4, 4.13241033e-4, 4.10672248e-4,
     4.08480950e-4, 4.07521131e-4, 4.06620556e-4, 4.05633162e-4, 4.04394743e-4, 4.02768826e-4, 4.01579752e-4, 4.00290096e-4, 3.98392706e-4, 3.96699174e-4,
     3.94234388e-4, 3.91851912e-4, 3.89154416e-4, 3.86231981e-4, 3.83777232e-4, 3.81168698e-4, 3.78433645e-4, 3.76298344e-4, 3.74215227e-4, 3.72197441e-4,
     3.70207284e-4, 3.68657365e-4, 3.67545261e-4, 3.66319403e-4, 3.65206679e-4, 3.63557693e-4, 3.62294538e-4, 3.61179418e-4, 3.58947889e-4, 3.56599608e-4,
     3.53923148e-4, 3.51224187e-4, 3.48070600e-4, 3.44364936e-4, 3.40781769e-4, 3.37471420e-4, 3.33912295e-4, 3.30235172e-4, 3.26561393e-4, 3.23282422e-4,
     3.20050433e-4, 3.17180477e-4, 3.14219227e-4, 3.11600742e-4, 3.09493381e-4, 3.07471257e-4, 3.05484183e-4, 3.03684277e-4, 3.01998980e-4, 3.00340107e-4,
     2.98343918e-4, 2.96867620e-4, 2.94945050e-4, 2.93219375e-4, 2.91688135e-4, 2.89922856e-4, 2.87800025e-4, 2.85461354e-4, 2.83757039e-4, 2.81377961e-4,
     2.79552629e-4, 2.77642310e-4, 2.75503391e-4, 2.73875335e-4, 2.72347155e-4, 2.70900748e-4, 2.69293738e-4, 2.67858074e-4, 2.66460300e-4, 2.65251110e-4,
     2.64275292e-4, 2.63248240e-4, 2.62327379e-4, 2.61579696e-4, 2.60818456e-4, 2.60002965e-4, 2.58705340e-4, 2.57752670e-4, 2.56913992e-4, 2.56361785e-4,
     2.55265269e-4, 2.54289030e-4, 2.53369783e-4, 2.52709129e-4, 2.51737103e-4, 2.50397300e-4, 2.48911620e-4, 2.47940106e-4, 2.46958831e-4, 2.46000402e-4,
     2.44657610e-4, 2.42959871e-4, 2.41646677e-4, 2.40071179e-4, 2.38572274e-4, 2.37125904e-4, 2.35921992e-4, 2.34980504e-4, 2.34028683e-4, 2.32572352e-4,
     2.31160035e-4, 2.29885626e-4, 2.28995291e-4, 2.28016467e-4, 2.26725040e-4, 2.25484801e-4, 2.23893466e-4, 2.22942513e-4, 2.21397987e-4, 2.19729430e-4,
     2.18184785e-4, 2.16237560e-4, 2.14867776e-4, 2.13119570e-4, 2.11095792e-4, 2.09431424e-4, 2.07369328e-4, 2.05197478e-4, 2.02992714e-4, 2.00995783e-4,
     1.99040628e-4, 1.97330365e-4, 1.95128937e-4, 1.92897050e-4, 1.90828812e-4, 1.88848992e-4, 1.86812770e-4, 1.85107403e-4, 1.83851878e-4, 1.81970227e-4,
     1.80471463e-4, 1.78981517e-4, 1.77286120e-4, 1.76392581e-4, 1.75098596e-4, 1.74137351e-4, 1.73190515e-4, 1.72166132e-4, 1.71231374e-4, 1.70366007e-4,
     1.69568760e-4, 1.68589267e-4, 1.67791363e-4, 1.66897578e-4, 1.65888439e-4, 1.65083935e-4, 1.64157795e-4, 1.63088458e-4, 1.62054638e-4, 1.60921204e-4,
     1.59597998e-4, 1.58591166e-4, 1.57448762e-4, 1.55975979e-4, 1.54945634e-4, 1.53768607e-4, 1.52490862e-4, 1.51237859e-4, 1.50099128e-4, 1.49111429e-4,
     1.47998498e-4, 1.47145525e-4, 1.46192385e-4, 1.45307609e-4, 1.44594903e-4, 1.43689016e-4, 1.42978840e-4, 1.42385947e-4, 1.41618642e-4, 1.41296459e-4,
     1.40736609e-4, 1.40128598e-4, 1.39496377e-4, 1.38800776e-4, 1.38248154e-4, 1.37619848e-4, 1.37121202e-4, 1.36531523e-4, 1.35719098e-4, 1.35125606e-4,
     1.34189288e-4, 1.33450058e-4, 1.32766920e-4, 1.31885017e-4, 1.31122573e-4, 1.30328488e-4, 1.29426836e-4, 1.28660144e-4, 1.27935543e-4, 1.27252022e-4,
     1.26480762e-4, 1.25668248e-4, 1.24922118e-4, 1.24410907e-4, 1.23808009e-4, 1.23163911e-4, 1.22517652e-4, 1.21996101e-4, 1.21426745e-4, 1.20795188e-4,
     1.20351946e-4, 1.19824156e-4, 1.19435754e-4, 1.19034938e-4, 1.18434743e-4, 1.18058836e-4, 1.17753346e-4, 1.17341202e-4, 1.17001634e-4, 1.16538566e-4,
     1.16315189e-4, 1.16025052e-4, 1.15877025e-4, 1.15561180e-4, 1.15133517e-4, 1.14836326e-4, 1.14582204e-4, 1.14243860e-4, 1.14197033e-4, 1.13911263e-4,
     1.13780203e-4, 1.13241434e-4, 1.12847757e-4, 1.12518590e-4, 1.12240191e-4, 1.12090644e-4, 1.11898539e-4, 1.11697770e-4, 1.11478439e-4, 1.11055915e-4,
     1.11038903e-4, 1.10892821e-4, 1.10844792e-4, 1.10832493e-4, 1.10685454e-4, 1.10585787e-4, 1.10523459e-4, 1.10538204e-4, 1.10352101e-4, 1.10342922e-4,
     1.10281513e-4, 1.10136790e-4, 1.10034253e-4, 1.10161106e-4, 1.10100572e-4, 1.10110098e-4, 1.10080592e-4, 1.09984173e-4, 1.09845052e-4, 1.09851854e-4,
     1.09905915e-4, 1.10028217e-4, 1.10035987e-4, 1.09960749e-4, 1.09766410e-4, 1.09785602e-4, 1.09799006e-4, 1.09935387e-4, 1.09996957e-4, 1.10015594e-4,
     1.09970966e-4, 1.09805155e-4, 1.09719944e-4, 1.09678392e-4, 1.09725399e-4, 1.09808313e-4, 1.09636200e-4, 1.09535133e-4, 1.09295215e-4, 1.09169349e-4,
     1.09191375e-4, 1.09233527e-4, 1.09205113e-4, 1.09120552e-4, 1.09119671e-4, 1.08973612e-4, 1.08924114e-4, 1.08891026e-4, 1.08907229e-4, 1.08986240e-4,
     1.08900531e-4, 1.08753082e-4, 1.08703896e-4, 1.08787039e-4, 1.08795093e-4, 1.08699656e-4, 1.08837754e-4, 1.08804051e-4, 1.08825203e-4, 1.08836245e-4,
     1.08750606e-4, 1.08911309e-4, 1.08993157e-4, 1.09019453e-4, 1.08906882e-4, 1.08830389e-4, 1.08786024e-4, 1.08648784e-4, 1.08543831e-4, 1.08526696e-4,
     1.08346842e-4, 1.08177325e-4, 1.07894351e-4, 1.07585013e-4, 1.07355722e-4, 1.07161797e-4, 1.07018269e-4, 1.06701081e-4, 1.06490135e-4, 1.06148921e-4,
     1.05780669e-4, 1.05519970e-4, 1.05237804e-4, 1.05018149e-4, 1.04960672e-4, 1.04759539e-4, 1.04540409e-4, 1.04335532e-4, 1.04248329e-4, 1.04080729e-4,
     1.04091425e-4, 1.04055083e-4, 1.04023531e-4, 1.03982955e-4, 1.03897792e-4, 1.03664565e-4, 1.03597003e-4, 1.03534202e-4, 1.03538706e-4, 1.03610086e-4,
     1.03562393e-4, 1.03382849e-4, 1.03208570e-4, 1.03080933e-4, 1.03024256e-4, 1.03051061e-4, 1.03046453e-4, 1.02968720e-4, 1.02813488e-4, 1.02603014e-4,
     1.02234419e-4, 1.02028022e-4, 1.01905582e-4, 1.01804529e-4, 1.01678913e-4, 1.01525915e-4, 1.01235735e-4, 1.01002510e-4, 1.00753430e-4, 1.00514078e-4,
     1.00437680e-4, 1.00426425e-4, 1.00342652e-4, 1.00324563e-4, 1.00193735e-4, 9.99614188e-5, 9.97959450e-5, 9.96635685e-5, 9.95680633e-5, 9.95295601e-5,
     9.95419478e-5, 9.93610409e-5, 9.92005015e-5, 9.89874238e-5, 9.87621992e-5, 9.86741797e-5, 9.86446806e-5, 9.87023188e-5, 9.86904988e-5, 9.87196887e-5,
     9.86048564e-5, 9.84391907e-5, 9.83676057e-5, 9.83207848e-5, 9.84386405e-5, 9.84955583e-5, 9.85247638e-5, 9.85229229e-5, 9.84406010e-5, 9.83803668e-5,
     9.83439478e-5, 9.83296520e-5, 9.84439321e-5, 9.85235870e-5, 9.85655663e-5, 9.84966496e-5, 9.84488197e-5, 9.84056559e-5, 9.83206327e-5, 9.81939782e-5,
     9.83002829e-5, 9.82246562e-5, 9.82700124e-5, 9.83010388e-5, 9.82759501e-5, 9.81680327e-5, 9.80868465e-5, 9.80081521e-5, 9.80084773e-5, 9.80201149e-5,
     9.81387400e-5, 9.80116911e-5, 9.80144351e-5, 9.79589533e-5, 9.78046061e-5, 9.76346492e-5, 9.76720556e-5, 9.76000615e-5, 9.76940237e-5, 9.76700649e-5,
     9.76436554e-5, 9.76217315e-5, 9.75644547e-5, 9.74582521e-5, 9.73723677e-5, 9.75158846e-5, 9.76063227e-5, 9.76576561e-5, 9.78378613e-5, 9.78764318e-5,
     9.78355668e-5, 9.79053556e-5, 9.78971200e-5, 9.79906854e-5, 9.80966670e-5, 9.83212609e-5, 9.83661316e-5, 9.85426549e-5, 9.86024674e-5, 9.86590075e-5,
     9.86276988e-5, 9.87354712e-5, 9.87754605e-5, 9.89135932e-5, 9.91850242e-5, 9.94027993e-5, 9.94630471e-5, 9.95331627e-5, 9.95927355e-5, 9.96449255e-5,
     9.96738155e-5, 9.97903713e-5, 9.99867011e-5, 1.00180641e-4, 1.00305263e-4, 1.00350675e-4, 1.00337790e-4, 1.00368296e-4, 1.00386643e-4, 1.00381686e-4,
     1.00519126e-4, 1.00708351e-4, 1.00900376e-4, 1.00985285e-4, 1.01155946e-4, 1.01232434e-4, 1.01246029e-4, 1.01294137e-4, 1.01423141e-4, 1.01638838e-4,
     1.01812099e-4, 1.01937894e-4, 1.02167422e-4, 1.02253729e-4, 1.02424793e-4, 1.02521930e-4, 1.02690121e-4, 1.02839311e-4, 1.03089113e-4, 1.03250475e-4,
     1.03431580e-4, 1.03716033e-4, 1.03938607e-4, 1.04117939e-4, 1.04384106e-4, 1.04529677e-4, 1.04654582e-4, 1.04814310e-4, 1.05092692e-4, 1.05413737e-4,
     1.05740355e-4, 1.06083661e-4, 1.06308365e-4, 1.06565438e-4, 1.06629335e-4, 1.06744210e-4, 1.06882442e-4, 1.07187881e-4, 1.07381175e-4, 1.07571985e-4,
     1.07752337e-4, 1.07958784e-4, 1.08010861e-4, 1.08075778e-4, 1.08124883e-4, 1.08257449e-4, 1.08419171e-4, 1.08593318e-4, 1.08694174e-4, 1.08884465e-4,
     1.09042040e-4, 1.09206302e-4, 1.09350437e-4, 1.09438760e-4, 1.09607358e-4, 1.09782787e-4, 1.09955758e-4, 1.10091358e-4, 1.10285940e-4, 1.10491229e-4,
     1.10742419e-4, 1.10907433e-4, 1.10989206e-4, 1.10994674e-4, 1.11224947e-4, 1.11275260e-4, 1.11379372e-4, 1.11598462e-4, 1.11899867e-4, 1.12232203e-4,
     1.12446731e-4, 1.12641100e-4, 1.12842223e-4, 1.13075513e-4, 1.13342586e-4, 1.13508149e-4, 1.13820379e-4, 1.14203203e-4, 1.14490663e-4, 1.14760307e-4,
     1.15092754e-4, 1.15469416e-4, 1.15670149e-4, 1.16009156e-4, 1.16271820e-4, 1.16580050e-4, 1.16902429e-4, 1.17308313e-4, 1.17666409e-4, 1.18021025e-4,
     1.18417002e-4, 1.18680029e-4, 1.18865319e-4, 1.19262047e-4, 1.19437496e-4, 1.19686615e-4, 1.19978762e-4, 1.20252739e-4, 1.20473473e-4, 1.20770056e-4,
     1.21090033e-4, 1.21331608e-4, 1.21505469e-4, 1.21699316e-4, 1.21782412e-4, 1.21974020e-4, 1.22103237e-4, 1.22231939e-4, 1.22420471e-4, 1.22581870e-4,
     1.22711136e-4, 1.22866496e-4, 1.22913279e-4, 1.23115711e-4, 1.23199889e-4, 1.23323819e-4, 1.23461205e-4, 1.23642315e-4, 1.23796507e-4, 1.24020468e-4,
     1.24238753e-4, 1.24422114e-4, 1.24625589e-4, 1.24860150e-4, 1.24989076e-4, 1.25232771e-4, 1.25461347e-4, 1.25751417e-4, 1.26039690e-4, 1.26363122e-4,
     1.26586473e-4, 1.26874300e-4, 1.27211480e-4, 1.27449080e-4, 1.27810195e-4, 1.28112243e-4, 1.28496922e-4, 1.28769318e-4, 1.28988874e-4, 1.29308334e-4,
     1.29681533e-4, 1.30081273e-4, 1.30484797e-4, 1.30834377e-4, 1.31195470e-4, 1.31542514e-4, 1.31971306e-4, 1.32327207e-4, 1.32671516e-4, 1.33163904e-4,
     1.33602896e-4, 1.34007504e-4, 1.34403545e-4, 1.34785880e-4, 1.35232311e-4, 1.35920055e-4, 1.36371634e-4, 1.36941354e-4, 1.37497873e-4, 1.38329979e-4,
     1.38981090e-4, 1.39719676e-4, 1.40596483e-4, 1.41573836e-4, 1.42701064e-4, 1.43912407e-4, 1.45127755e-4, 1.46710702e-4, 1.48201619e-4, 1.49963009e-4,
     1.51809861e-4, 1.53750193e-4, 1.55815317e-4, 1.57693897e-4, 1.59629003e-4, 1.61541153e-4, 1.63306339e-4, 1.64888939e-4, 1.66286077e-4, 1.67668695e-4,
     1.68815611e-4, 1.69629768e-4, 1.70307707e-4, 1.70592642e-4, 1.70837390e-4, 1.70905017e-4, 1.70730284e-4, 1.70560321e-4, 1.70279811e-4, 1.69980969e-4,
     1.69683568e-4, 1.69203059e-4, 1.68858534e-4, 1.68509619e-4, 1.68194840e-4, 1.67681333e-4, 1.67199028e-4, 1.66717819e-4, 1.66199767e-4, 1.65735840e-4,
     1.65420369e-4, 1.65032314e-4, 1.64890865e-4, 1.64715724e-4, 1.64497541e-4, 1.64412240e-4, 1.64324793e-4, 1.64276994e-4, 1.64155537e-4, 1.64020542e-4,
     1.64016963e-4, 1.63898890e-4, 1.64032362e-4, 1.64045267e-4, 1.64195414e-4, 1.64285742e-4, 1.64460407e-4, 1.64661763e-4, 1.64738587e-4, 1.64862909e-4,
     1.65076816e-4, 1.65024997e-4, 1.65162453e-4, 1.65237474e-4, 1.65562552e-4, 1.65709225e-4, 1.66142295e-4, 1.66413982e-4, 1.66803643e-4, 1.67161297e-4,
     1.67479188e-4, 1.67667097e-4, 1.67864791e-4, 1.68139958e-4, 1.68233347e-4, 1.68368232e-4, 1.68734109e-4, 1.68895921e-4, 1.69306212e-4, 1.69651966e-4,
     1.69946857e-4, 1.70444809e-4, 1.70732732e-4, 1.71135157e-4, 1.71231516e-4, 1.71550211e-4, 1.71696417e-4, 1.71827084e-4, 1.72142950e-4, 1.72525391e-4,
     1.72856514e-4, 1.73390715e-4, 1.73767822e-4, 1.74328441e-4, 1.74893554e-4, 1.75395641e-4, 1.75742115e-4, 1.76164298e-4, 1.76419331e-4, 1.76746627e-4,
     1.76913681e-4, 1.77372254e-4, 1.77897899e-4, 1.78453513e-4, 1.79169646e-4, 1.79727818e-4, 1.80409802e-4, 1.81135436e-4, 1.81693575e-4, 1.82304987e-4,
     1.82840211e-4, 1.83218701e-4, 1.83686875e-4, 1.84158245e-4, 1.84720710e-4, 1.85261670e-4, 1.86044882e-4, 1.86910627e-4, 1.87744111e-4, 1.88640197e-4,
     1.89536957e-4, 1.90373887e-4, 1.91044376e-4, 1.91770554e-4, 1.92334227e-4, 1.93008156e-4, 1.93752580e-4, 1.94471581e-4, 1.95272205e-4, 1.96159946e-4,
     1.97179878e-4, 1.98147132e-4, 1.99205315e-4, 2.00512171e-4, 2.01401775e-4, 2.02350741e-4, 2.03121965e-4, 2.03880283e-4, 2.04553652e-4, 2.05224749e-4,
     2.05896678e-4, 2.06741041e-4, 2.07620395e-4, 2.08540576e-4, 2.09373829e-4, 2.10385794e-4, 2.11266605e-4, 2.12109107e-4, 2.12682109e-4, 2.13305338e-4,
     2.13902214e-4, 2.14484164e-4, 2.14843762e-4, 2.15242449e-4, 2.15718844e-4, 2.16421684e-4, 2.17195290e-4, 2.18005590e-4, 2.18962077e-4, 2.19991471e-4,
     2.20919919e-4, 2.21546974e-4, 2.22053928e-4, 2.22555080e-4, 2.23237798e-4, 2.23532485e-4, 2.23996613e-4, 2.24501430e-4, 2.25173147e-4, 2.25737767e-4,
     2.26423767e-4, 2.27365037e-4, 2.28449699e-4, 2.29400621e-4, 2.30377333e-4, 2.31057918e-4, 2.31964152e-4, 2.32623742e-4, 2.33166196e-4, 2.33650080e-4,
     2.34105936e-4, 2.34878659e-4, 2.35394209e-4, 2.36099367e-4, 2.37123530e-4, 2.38223666e-4, 2.39309613e-4, 2.40405595e-4, 2.41495657e-4, 2.42458496e-4,
     2.43313474e-4, 2.44175726e-4, 2.44734401e-4, 2.45486851e-4, 2.45993980e-4, 2.46800743e-4, 2.47463201e-4, 2.48445219e-4, 2.49488279e-4, 2.50569021e-4,
     2.51747840e-4, 2.52951109e-4, 2.54188859e-4, 2.55386833e-4, 2.56298369e-4, 2.57479006e-4, 2.58213101e-4, 2.59065147e-4, 2.60018887e-4, 2.61076885e-4,
     2.62283407e-4, 2.63904910e-4, 2.65792030e-4, 2.67956880e-4, 2.70494493e-4, 2.73225853e-4, 2.76170072e-4, 2.79055057e-4, 2.81984548e-4, 2.88500954e-4,
     2.89109910e-4, 2.93129400e-4, 2.92536512e-4, 2.92536512e-4] #[uW m-2 xps-1] 10 sensitivities per line, 1024 total


def wavelengthSpectrumAnddownwellingSpectralFlux(fileLocation):
    '''
    This function will format the source JSON file including multiple JSON objects
    into a file of JSON array
    '''
    with open(fileLocation, 'r') as fileHandler:
        tempList, wvl_lgr, spectrum, k, writeToWavelength =\
            fileHandler.read().split('\n'), [], [[]], 0, True

        for i in range(len(tempList)):
            if "wavelength" in tempList[i] and writeToWavelength:
                wvl_lgr.append(
                    float(tempList[i][tempList[i].find(':') + 1: -2]) * 1e-9)
            if "wavelength" not in tempList[i] and "wavelength" in tempList[i - 4]\
                    and "band" not in tempList[i] and "," not in tempList[i]:
                writeToWavelength = False
                spectrum.append([])
                k += 1
            if "spectrum" in tempList[i]:
                spectrum[k].append(
                    float(tempList[i][tempList[i].find(':') + 1: -2]))

        spectrum.remove([])

    wvl_ntf  = [np.average([wvl_lgr[i], wvl_lgr[i+1]]) for i in range(len(wvl_lgr)-1)]
    delta = [wvl_ntf[i+1] - wvl_ntf[i] for i in range(len(wvl_ntf) - 1)]
    delta.insert(0, 2*(wvl_ntf[0] - wvl_lgr[0]))
    delta.insert(-1, 2*(wvl_lgr[-1] - wvl_ntf[-1]))

    Spectrometer_Integration_Time_In_Microseconds = 5000.0
    Spectrometer_Integration_Time = Spectrometer_Integration_Time_In_Microseconds / 1.0e6

    Area= np.pi * (3900.0*1.0e-6) ** 2 / 4.0

    #General formula used in calculating downwelling spectral flux:
    #Downwelling Spectral Flux = flx_sns * spectrum / bandwidth / area / time
    downwellingSpectralFlux = np.array(_FLX_SNS) * np.array(spectrum) * 1.0e-6 / np.array(delta) / Area / Spectrometer_Integration_Time

    return wvl_lgr, spectrum, downwellingSpectralFlux


def JSONHandler(fileLocation):
    '''
    Main JSON handler, write JSON file to a Python list with standard JSON module
    '''
    wavelength, spectrum, downwellingSpectralFlux = wavelengthSpectrumAnddownwellingSpectralFlux(
        fileLocation)
    with open(fileLocation, 'r') as fileHandler:
        return json.loads(fileHandler.read()), wavelength, spectrum, downwellingSpectralFlux


def renameTheValue(name):
    '''
    Rename the value so it becomes legal in netCDF
    '''
    if type(name) is unicode:
        name = name.encode('ascii', 'ignore')
    if name in _UNIT_DICTIONARY:
        name = _UNIT_DICTIONARY[name]
    elif name in _NAMES:
        name = _NAMES[name]

    return name.replace(" ", "_")


def getSpectrometerInformation(arrayOfJSON):
    '''
    Collect information from spectrometer with special care
    '''
    maxFixedIntensity = [int(intensityMembers["spectrometer"]["maxFixedIntensity"]) for intensityMembers in
                         arrayOfJSON]
    integrationTime = [int(integrateMembers["spectrometer"]["integration time in ?s"]) for integrateMembers in
                       arrayOfJSON]

    return maxFixedIntensity, integrationTime


def getListOfValue(arrayOfJSON, dataName):
    '''
    Collect data from JSON objects which have "value" member
    '''
    return [float(valueMembers[dataName]['value'].encode('ascii', 'ignore')) for valueMembers in arrayOfJSON]


def getListOfRawValue(arrayOfJSON, dataName):
    '''
    Collect data from JSON objects which have "rawValue" member
    '''
    return [float(valueMembers[dataName]['rawValue'].encode('ascii', 'ignore')) for valueMembers in arrayOfJSON]


def _timeStamp():
    '''
    Record the time the script is triggered
    '''
    return time.strftime("%a %b %d %H:%M:%S %Y",  time.localtime(int(time.time())))


def translateTime(timeString):
    '''
    Translate the time the metadata included as the days offset to the basetime.
    '''
    timeUnpack = datetime.strptime(timeString, "%Y.%m.%d-%H:%M:%S").timetuple()
    timeSplit = date(year=timeUnpack.tm_year, month=timeUnpack.tm_mon,
                     day=timeUnpack.tm_mday) - _UNIX_BASETIME

    return (timeSplit.total_seconds() + timeUnpack.tm_hour * 3600.0 + timeUnpack.tm_min * 60.0 +
            timeUnpack.tm_sec) / (3600.0 * 24.0)


def main(JSONArray, outputFileName, wavelength=None, spectrum=None, downwellingSpectralFlux=None, recordTime=None, commandLine=None):
    '''
    Main netCDF handler, write data to the netCDF file indicated.
    '''
    netCDFHandler         = Dataset(outputFileName, 'w', format='NETCDF4')
    dataMemberList        = [JSONMembers[u"environment_sensor_set_reading"]
                      for JSONMembers in JSONArray]
    timeList              = [translateTime(JSONMembers[u'timestamp'])
                for JSONMembers in dataMemberList]
    timeStampList         = [JSONMembers[u'timestamp']
                     for JSONMembers in dataMemberList]
    timeDimension         = netCDFHandler.createDimension("time", None)
    tempTimeStampVariable = netCDFHandler.createVariable("timestamp", str, ("time",), chunksizes=(1,))
    for i in range(len(timeStampList)):
        tempTimeStampVariable[i] = timeStampList[i]
    tempTimeVariable      = netCDFHandler.createVariable(
        'time', 'f8', ('time',), chunksizes=(1,))
    for i in range(len(timeList)):  # Assign Times
        tempTimeVariable[i] = timeList[i]
    setattr(tempTimeVariable, "units",    "days since 1970-01-01 00:00:00")
    setattr(tempTimeVariable, "calender", "gregorian")

    for data in dataMemberList[0]:
        if data != 'spectrometer' and type(dataMemberList[0][data]) not in (str, unicode):
            tempVariable = netCDFHandler.createVariable(
                renameTheValue(data), 'f4', ('time',))
            tempVariable[:] = getListOfValue(
                dataMemberList, data)  # Assign "values"
            if 'unit' in dataMemberList[0][data]:  # Assign Units
                setattr(tempVariable, 'units', _UNIT_DICTIONARY[
                        dataMemberList[0][data]['unit']])
            if 'rawValue' in dataMemberList[0][data]:  # Assign "rawValues"
                netCDFHandler.createVariable(renameTheValue(data) + '_rawValue', 'f4', ('time',))[:] =\
                    getListOfRawValue(dataMemberList, data)
        elif type(dataMemberList[0][data]) in (str, unicode) and data != "timestamp":
            netCDFHandler.createVariable(renameTheValue(data), str)[
                0] = dataMemberList[0][data]

        if data == 'spectrometer':  # Special care for spectrometers :)
            netCDFHandler.createVariable('Spectrometer_maxFixedIntensity', 'f4', ('time',))[:] =\
                getSpectrometerInformation(dataMemberList)[0]
            netCDFHandler.createVariable('Spectrometer_Integration_Time_In_Microseconds', 'f4', ('time',))[:] =\
                getSpectrometerInformation(dataMemberList)[1]

    if wavelength and spectrum:
        netCDFHandler.createDimension("wvl_lgr", len(wavelength))
        netCDFHandler.createVariable("wvl_lgr", 'f4', ('wvl_lgr',))[
            :] = wavelength
        setattr(netCDFHandler.variables['wvl_lgr'], 'units', 'meter')
        setattr(netCDFHandler.variables['wvl_lgr'], 'long_name', 'Wavelength')
        netCDFHandler.createVariable("spectrum", 'f4', ('time', 'wvl_lgr'))[
            :, :] = spectrum

    wvl_ntf  = [np.average([wavelength[i], wavelength[i+1]]) for i in range(len(wavelength)-1)]
    delta = [wvl_ntf[i+1] - wvl_ntf[i] for i in range(len(wvl_ntf) - 1)]
    delta.insert(0, 2*(wvl_ntf[0] - wavelength[0]))
    delta.insert(-1, 2*(wavelength[-1] - wvl_ntf[-1]))

    # Add data from terraref.nco
    netCDFHandler.createVariable("wvl_dlt", 'f8', ("wvl_lgr",))[:] = delta
    setattr(netCDFHandler.variables['wvl_dlt'], 'units', 'meter')
    setattr(netCDFHandler.variables['wvl_dlt'], 'notes',"Best information available is that bandwidth is uniform 0.4 nm across all channels")
    setattr(netCDFHandler.variables['wvl_dlt'], 'long_name', "Bandwidth of environmental sensor")

    netCDFHandler.createVariable("flx_sns", "f4", ("wvl_lgr",))[:] = [members/1.0e-6 for members in _FLX_SNS]
    setattr(netCDFHandler.variables['flx_sns'],'units', 'watt meter-2 count-1')
    setattr(netCDFHandler.variables['flx_sns'],'long_name','Flux sensitivity of each band (irradiance per count)')
    setattr(netCDFHandler.variables['flx_sns'], 'provenance', "EnvironmentalLogger calibration information from file S05673_08062015.IrradCal provided by TinoDornbusch and discussed here: https://github.com/terraref/reference-data/issues/30#issuecomment-217518434")

    netCDFHandler.createVariable("flx_spc_dwn", 'f4', ('time','wvl_lgr'))[:,:] = downwellingSpectralFlux
    setattr(netCDFHandler.variables['flx_spc_dwn'],'units', 'watt meter-2 meter-1')
    setattr(netCDFHandler.variables['flx_spc_dwn'], 'long_name', 'Downwelling Spectral Irradiance')

    # Downwelling Flux = summation of (delta lambda(_wvl_dlt) * downwellingSpectralFlux)
    netCDFHandler.createVariable("flx_dwn", 'f4').assignValue(np.sum(downwellingSpectralFlux))
    setattr(netCDFHandler.variables["flx_dwn"], "units", "watt meter-2")
    setattr(netCDFHandler.variables['flx_dwn'], 'long_name', 'Downwelling Irradiance')

    #Other Constants used in calculation
    #Integration Time
    netCDFHandler.createVariable("int_tme", 'f4').assignValue(5000.0/1.0e-6)
    setattr(netCDFHandler.variables["flx_dwn"], "units", "second")
    setattr(netCDFHandler.variables['flx_dwn'], 'long_name', 'Spectraometer integration time')

    #Spectrometer Area
    netCDFHandler.createVariable("area", "f4").assignValue(np.pi * (3900.0*1.0e-6) ** 2 / 4.0)
    setattr(netCDFHandler.variables["area"], "units", "meter2")
    setattr(netCDFHandler.variables['area'], 'long_name', 'Spectraometer Area')



    netCDFHandler.history = recordTime + ': python ' + commandLine
    netCDFHandler.close()

if __name__ == '__main__':
    fileInputLocation, fileOutputLocation = sys.argv[1], sys.argv[2]
    if not os.path.exists(fileOutputLocation) and not fileOutputLocation.endswith('.nc'):
        os.mkdir(fileOutputLocation)  # Create folder

    if not os.path.isdir(fileInputLocation) or fileOutputLocation.endswith('.nc'):
        print "Processing", fileInputLocation + '....'
        tempJSONMasterList, wavelength, spectrum, downwellingSpectralFlux = JSONHandler(
            fileInputLocation)
        if not os.path.isdir(fileOutputLocation):
            main(tempJSONMasterList, fileOutputLocation, wavelength, spectrum, downwellingSpectralFlux,
                 _timeStamp(), sys.argv[1] + ' ' + sys.argv[2])
        else:
            outputFileName = os.path.split(fileInputLocation)[-1]
            main(tempJSONMasterList, os.path.join(fileOutputLocation,
                                                  outputFileName.strip('.json') + '.nc'), wavelength, spectrum, downwellingSpectralFlux,
                 _timeStamp(), sys.argv[1] + ' ' + sys.argv[2])
    else:  # Read and Export netCDF to folder
        for filePath, fileDirectory, fileName in os.walk(fileInputLocation):
            for members in fileName:
                if os.path.join(filePath, members).endswith('.json'):
                    print "Processing", members + '....'
                    outputFileName = members.strip('.json') + '.nc'
                    tempJSONMasterList, wavelength, spectrum, downwellingSpectralFlux = JSONHandler(
                        os.path.join(filePath, members))
                    main(tempJSONMasterList, os.path.join(
                        fileOutputLocation, outputFileName),
                        wavelength, spectrum, downwellingSpectralFlux, _timeStamp(), sys.argv[1] + ' ' + sys.argv[2])
