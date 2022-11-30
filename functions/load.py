def findfiles(directory,fileformat,daylist):
	# find all files with fileformat under directory
	# if daylist is specified, only find files in certain days
	# ex) daylist = [3,5,6]: search files in directory/Day3, directory/Day5, directory/Day6
	# if daylist is empty, search all files in directory
	import os
	def findday(file):
		if 'Day' in file:
			import re
			return int(re.split('Day|_',os.path.basename(os.path.dirname(file)))[1])
		else:
			return 0

	files = [os.path.join(root, name)
			 for root, dirs, files in os.walk(directory)
			 for name in files if os.path.splitext(name)[1] in fileformat]
	files = sorted(files, key=findday)
	days = [findday(f) for f in files]
	files = [x for x,y in zip(files,days) if y>0]
	days = [x for x in days if x>0]

	if len(daylist)>0:
		indaylist = [i for i,v in enumerate(days) if v in daylist]
		files = [files[i] for i in indaylist]
		days = [days[i] for i in indaylist]

	return files, days


def load_mat(filename):
	import numpy as np
	from scipy.io import loadmat, matlab

	"""
    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

	def _check_vars(d):
		"""
        Checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
		for key in d:
			if isinstance(d[key], matlab.mio5_params.mat_struct):
				d[key] = _todict(d[key])
			elif isinstance(d[key], np.ndarray):
				d[key] = _toarray(d[key])
		return d

	def _todict(matobj):
		"""
        A recursive function which constructs from matobjects nested dictionaries
        """
		d = {}
		for strg in matobj._fieldnames:
			elem = matobj.__dict__[strg]
			if isinstance(elem, matlab.mio5_params.mat_struct):
				d[strg] = _todict(elem)
			elif isinstance(elem, np.ndarray):\
				d[strg] = _toarray(elem)
			else:
				d[strg] = elem
		return d
	def _toarray(ndarray):
		"""
        A recursive function which constructs ndarray from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
		if ndarray.dtype != 'float64':
			elem_list = []
			for sub_elem in ndarray:
				if isinstance(sub_elem, matlab.mio5_params.mat_struct):
					elem_list.append(_todict(sub_elem))
				elif isinstance(sub_elem, np.ndarray):
					elem_list.append(_toarray(sub_elem))
				else:
					elem_list.append(sub_elem)
			return np.array(elem_list)
		else:
			return ndarray

	data = loadmat(filename, struct_as_record=False, squeeze_me=True)
	return _check_vars(data)

def load_doric(filename,version,datanames,datanames_new):
	# loading doric file directly, this is modified version of code from doric
	#
	# filename: doric file (.doric) name
	# version: doric neuroscience studio version; this needs to be 5 for now;
	# this is for generalizing function across files that recorded using different version of software
	# datanames: how each data is saved in raw doric file
	# datanames_new: how we want to call them - for generalization with photometry file from pyphotometry system
	import h5py
	import numpy as np

	def ish5dataset(item):
		return isinstance(item, h5py.Dataset)

	def h5getDatasetR(item, leading=''):
		r = []
		for key in item:
			# First have to check if the next layer is a dataset or not
			firstkey = list(item[key].keys())[0]
			if ish5dataset(item[key][firstkey]):
				r = r + [{'Name': leading + '_' + key, 'Data':
					[{'Name': k, 'Data': np.array(item[key][k]),
					  'DataInfo': {atrib: item[key][k].attrs[atrib] for atrib in item[key][k].attrs}} for k in
					 item[key]]}]
			else:
				r = r + h5getDatasetR(item[key], leading + '_' + key)

		return r

	# Extact Data from a doric file
	def ExtractDataAcquisition(filename, version):
		output = {}
		with h5py.File(filename, 'r') as h:
			# print(filename)
			if version == 5:
				return h5getDatasetR(h['Traces'], filename)
			elif version == 6:
				return h5getDatasetR(h['DataAcquisition'], filename)

	doricfile = [data["Data"][0] for data in ExtractDataAcquisition(filename,version)]
	data = {}
	if version == 5:
		data['time'] = [data["Data"] for data in doricfile if data["Name"] == 'Console_time(s)'][0]
	for i, v in enumerate(datanames):
		if version == 5:
			data[datanames_new[i]] = [data["Data"] for data in doricfile if data["Name"] == v][0]

	return data

def load_ppd(filename,datanames,datanames_new):
	# loading pyphotometry file directly, this is modified version of code from pyphotometry
	#
	# filename: pyphotometry (.py) file name
	# datanames: how each data is saved in raw py file
	# datanames_new: how we want to call them - for generalization with photometry file from doric system

	import json
	import numpy as np
	from scipy.signal import butter, filtfilt

	def import_ppd(file_path, low_pass=20, high_pass=0.01):
		'''Function to import pyPhotometry binary data files into Python. The high_pass
        and low_pass arguments determine the frequency in Hz of highpass and lowpass
        filtering applied to the filtered analog signals. To disable highpass or lowpass
        filtering set the respective argument to None.  Returns a dictionary with the
        following items:
            'subject_ID'    - Subject ID
            'date_time'     - Recording start date and time (ISO 8601 format string)
            'mode'          - Acquisition mode
            'sampling_rate' - Sampling rate (Hz)
            'LED_current'   - Current for LEDs 1 and 2 (mA)
            'version'       - Version number of pyPhotometry
            'analog_1'      - Raw analog signal 1 (volts)
            'analog_2'      - Raw analog signal 2 (volts)
            'analog_1_filt' - Filtered analog signal 1 (volts)
            'analog_2_filt' - Filtered analog signal 2 (volts)
            'digital_1'     - Digital signal 1
            'digital_2'     - Digital signal 2
            'pulse_inds_1'  - Locations of rising edges on digital input 1 (samples).
            'pulse_inds_2'  - Locations of rising edges on digital input 2 (samples).
            'pulse_times_1' - Times of rising edges on digital input 1 (ms).
            'pulse_times_2' - Times of rising edges on digital input 2 (ms).
            'time'          - Time of each sample relative to start of recording (ms)
        '''
		with open(file_path, 'rb') as f:
			header_size = int.from_bytes(f.read(2), 'little')
			data_header = f.read(header_size)
			data = np.frombuffer(f.read(), dtype=np.dtype('<u2'))
		# Extract header information
		header_dict = json.loads(data_header)
		volts_per_division = header_dict['volts_per_division']
		sampling_rate = header_dict['sampling_rate']
		# Extract signals.
		analog = data >> 1  # Analog signal is most significant 15 bits.
		digital = ((data & 1) == 1).astype(int)  # Digital signal is least significant bit.
		# Alternating samples are signals 1 and 2.
		analog_1 = analog[::2] * volts_per_division[0]
		analog_2 = analog[1::2] * volts_per_division[1]
		digital_1 = digital[::2]
		digital_2 = digital[1::2]
		time = np.arange(analog_1.shape[0]) * 1000 / sampling_rate  # Time relative to start of recording (ms).
		# Filter signals with specified high and low pass frequencies (Hz).
		if low_pass and high_pass:
			b, a = butter(2, np.array([high_pass, low_pass]) / (0.5 * sampling_rate), 'bandpass')
		elif low_pass:
			b, a = butter(2, low_pass / (0.5 * sampling_rate), 'low')
		elif high_pass:
			b, a = butter(2, high_pass / (0.5 * sampling_rate), 'high')
		if low_pass or high_pass:
			analog_1_filt = filtfilt(b, a, analog_1)
			analog_2_filt = filtfilt(b, a, analog_2)
		else:
			analog_1_filt = analog_2_filt = None
		# Extract rising edges for digital inputs.
		pulse_inds_1 = 1 + np.where(np.diff(digital_1) == 1)[0]
		pulse_inds_2 = 1 + np.where(np.diff(digital_2) == 1)[0]
		pulse_times_1 = pulse_inds_1 * 1000 / sampling_rate
		pulse_times_2 = pulse_inds_2 * 1000 / sampling_rate
		# Return signals + header information as a dictionary.
		data_dict = {'analog_1': analog_1,
					 'analog_2': analog_2,
					 'analog_1_filt': analog_1_filt,
					 'analog_2_filt': analog_2_filt,
					 'digital_1': digital_1,
					 'digital_2': digital_2,
					 'pulse_inds_1': pulse_inds_1,
					 'pulse_inds_2': pulse_inds_2,
					 'pulse_times_1': pulse_times_1,
					 'pulse_times_2': pulse_times_2,
					 'time': time}
		data_dict.update(header_dict)
		return data_dict
	data_temp = import_ppd(filename)
	data = {}
	data['time'] = data_temp['time']
	for i,v in enumerate(datanames):
		data[datanames_new[i]] = data_temp[v]
	return data

def load_pickle(filename):
	# load pickle file
	import pickle
	objects = []
	with (open(filename, "rb")) as openfile:
		while True:
			try:
				objects.append(pickle.load(openfile))
			except EOFError:
				break
	return objects

def load_dandi_url(dandiset_id,animalname,daylist=None):
	# load url in dandi server for specific sessions of given animal
	# dandiset_id: dandi set id, which is given by dandi when you upload the data
	# animalname: name of animal, it needs to be the same with one on dandi server
	# daylist: indicator for session (the session name needs to contain 'DayXX')
	# if you don't have daylist, it will give urls for all sessions of given animal
	from dandi.dandiapi import DandiAPIClient
	import numpy as np
	import re

	url = []
	path = []
	with DandiAPIClient.for_dandi_instance("dandi") as client:
		dandiset = client.get_dandiset(dandiset_id, 'draft')
		for asset in dandiset.get_assets_by_glob(animalname):
			url = np.append(url, asset.get_content_url(follow_redirects=1, strip_query=True))
			path = np.append(path,asset.get_metadata().path)

	day = [int(re.split('.nwb|-',x.split('Day')[1])[0]) for x in path]
	if not daylist==None:
		url = [y for x,y in zip(day,url) if x in daylist]
		path = [y for x,y in zip(day,path) if x in daylist]
		day = [x for x in day if x in daylist]

	url = [y for x, y in sorted(zip(day, url))] # this is url of each session
	path = [y for x, y in sorted(zip(day, path))] # this is file name, which follows this: sub_(animalname)_ses_(sessionname)

	return url, path


def load_nwb(url,namespacepath,varlist):
	# load nwb file using url - this needs to be developed more to be able to load saved nwb file
	# url: dandi server url
	# namespacepath: path for namespace of namboodirilab extension file
	# varlist: the list of variable you want to get
	# (a,XX): variable XX from acquisition field
	# (p,XX): variable XX from processing field

	from pynwb import NWBHDF5IO,load_namespaces
	import numpy as np

	for ipath in namespacepath:
		load_namespaces(ipath)

	io = NWBHDF5IO(url, mode='r', driver='ros3')
	nwbfile = io.read()

	results = {}
	for i in varlist:
		if i[0] =='a': #acquisition
			fields = nwbfile.acquisition[i[1]]._get_fields()
			subnwb = nwbfile.acquisition[i[1]]
		elif i[0] == 'p': #processing
			fields = nwbfile.processing[i[1]][i[2]]._get_fields()
			subnwb = nwbfile.processing[i[1]][i[2]]
		else:
			fields = nwbfile.i[1]._get_fields()
			subnwb = nwbfile.i[1]
		results[i[-1]] = {}
		for f in fields:
			if not subnwb.fields.get(f)==None:
				results[i[-1]][f] = subnwb.fields.get(f)
				if not np.shape(results[i[-1]][f])==():
					results[i[-1]][f] = results[i[-1]][f][:]
	return results, nwbfile

