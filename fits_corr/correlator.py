# Imports
from libs import utils
from os import path, makedirs
from sys import stdout
from multiprocessing import Process, Queue, cpu_count
from skimage.measure import compare_ssim
import warnings
import numpy as np
import json
from collections import Iterable

class Correlator(object):
    
    def __init__(self):
        self.props = {
            'align_window': 250,
            'prmi_radius': 3,
            'histogram_bins': 256,
            'multi_process': False,
            'data_dir': None,
            'data_pair': None,
            'out_dir': './out/',
            'coalesce': False
        }
        self.keys = [
            ('prmi', self._prmi), 
            ('nmi', self._mi), 
            ('pcc', self._pcc), 
            ('ssim', self._ssim)
        ]
        
    def _config_data_pair_helper(self, value):
        """
        Internal helper method to validate value given by .json config
        file for single on and off pair task.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        if value:
            value_iter = isinstance(v, Iterable) and 0 < len(v) < 3
            assert value_iter, 'Error: Data Pair should be None or 2-element Array!'
            self.with_on_off_pair(*value)
        
        return self
            
        
    def load_config(self, config_path):
        """
        Load configuration from a JSON file specified by <config_path>.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        err_msg = 'Error: Specified Config File Doesn\'t Exist or not a JSON file.'
        assert path.isfile(config_path) and config_path.endswith('.json'), err_msg
        
        config = {
            'align_window': lambda v: self.set_align_window(v),
            'prmi_radius': lambda v: self.set_prmi_radius(v),
            'histogram_bins': lambda v: self.set_histogram_bins(v),
            'multi_process': lambda v: v and self.enable_multiprocessing(),
            'data_dir': lambda v: v and self.with_data_dir(v),
            'data_pair': self._config_data_pair_helper,
            'out_dir': lambda v: v and self.with_output_dir(v),
            'coalesce': lambda v: v and self.enable_coalescing()
        }
        
        with open(config_path, 'r') as fp:
            for key, value in json.load(fp).items():
                if key in config:
                    config[key](value)
        
        print '[INFO] Configuration Loaded Successfully!'
        
        return self
    
    def save_config(self, save_dir):
        """
        Saves the current configuration to the directory specified by
        <save_dir> into a file named _CONFIG.json.
        
        Input
        -----
        save_dir  : Save directory path.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        if not path.exists(save_dir):
            makedirs(save_dir)
        
        self._json_writer(self.props, path.join(save_dir, '_CONFIG.json'))
        
        return self
    
    def set_align_window(self, span):
        """
        Set the default search span used during alignment.
        
        Input
        -----
        span      : Alignment window width.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        self.props['align_window'] = span
        
        return self
    
    def set_prmi_radius(self, radius):
        """
        Set the neighborhood radius for PCA Regional Mutual Information.
        
        Input
        -----
        radius    : Neighborhood radius.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        self.props['prmi_radius'] = radius
        
        return self
    
    def set_histogram_bins(self, bins):
        """
        Set the total <bins> to use in building the joint histogram for
        mutual information.
        
        Input
        -----
        bins    : Total bins for each dimension. Type: int or 2-tuple. 
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        bins_int = bins.__class__.__name__ == 'int'
        bins_iter = (isinstance(bins, Iterable) and 0 < len(bins) < 3)
        
        assert bins_int or bins_iter, 'Error: Invalid Specification for Bins!'
        
        self.props['histogram_bins'] = bins
        
        return self
        
    def enable_multiprocessing(self):
        """
        Enable Multi-Threading during processing the .fits files. Only
        available if the correlator is assigned a directory as its
        primary task.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        if self.props['data_dir'] and not self.props['data_pair']:
            self.props['multi_process'] = True
        else:
            message = 'Multiprocessing not enabled. ' + \
                      'Check that data directory is assigned and data pair is unassigned.'
            warnings.warn(message, RuntimeWarning)
        
        return self
    
    def with_data_dir(self, dirname):
        """
        Assign directory of .fits files given by <dirname> as the primary 
        task of the correlator.
        
        Inputs
        ------
        dirname   : Path to the directory of .fits files.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        assert path.exists(dirname), 'Error: Specified Task Directory Not Found!'
        
        self.props['data_dir'] = dirname
        self.props['data_pair'] = None
        
        return self
        
    def with_on_off_pair(self, on_fp, off_fp):
        """
        Assign ON and OFF pair of .fits file as the primary task of the
        correlator. This takes priority over any directory assigned as
        the primary task.
        
        Inputs
        ------
        on_fp     : Path to the ON .fits file.
        off_fp    : Path to the OFF .fits file.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        assert path.isfile(on_fp), 'Error: ON File Doesn\'t Exist!'
        assert path.isfile(off_fp), 'Error: OFF File Doesn\'t Exist!'
        
        self.props['data_pair'] = (on_fp, off_fp)
        self.props['data_dir'] = None
        
        return self
    
    def with_output_dir(self, dirname):
        """
        Assign an output directory given by <dirname> to the correlator for
        dumping the .json files.
        
        Input
        -----
        dirname   : Path to the output directory. Will be created when run() is 
                    called if non-existent.
        
        Output
        ------
        Returns the current correlator (self).
        
        """        
        self.props['out_dir'] = dirname
        
        return self
    
    def enable_coalescing(self):
        """
        Coalesce output scores to a single .json file.
        
        Output
        ------
        Returns the current correlator (self).
        
        """
        self.props['coalesce'] = True
        
        if self.props['data_pair']:
            message = 'Coalescing is a non-effect with data pairs.'
            warnings.warn(message, RuntimeWarning)
        
        return self
    
    def _translate(self, on_data, off_data, shift):
        """
        Internal method to translate <on_data> against the <off_data> by
        amount given by <shift>.
        
        Inputs
        ------
        on_data   : ON data ndarray.
        off_data  : OFF data ndarray.
        shift     : Amount to translate. Can be negative.
        
        Outputs
        -------
        The shifted <on_data> and <off_data> pair.
        
        """
        if shift < 0:
            return (on_data[:, -shift:], off_data[:, :shift])
        elif shift > 0:
            return (on_data[:, :-shift], off_data[:, shift:])
        else:
            return (on_data, off_data)
    
    def _align(self, on_data, off_data):
        """
        Deprecated Method. Use _align_fast instead.
        """
        valid_window = 0 < self.props['align_window'] < on_data.shape[1]
        assert valid_window, 'Error: Invalid Window Size!'
        
        score = lambda a, b: np.sum((a-a.mean()) * (b-b.mean()))/(a.std() * b.std())
        sim = lambda z: ((1+z)/2).max()
        area = lambda w: on_data.shape[0] * (on_data.shape[1] - abs(w))
        
        lower = -int(self.props['align_window']/2)
        upper = 1-lower
        
        scores = np.vectorize(
            lambda w: score(*self._translate(on_data, off_data, w))/area(w)
        )(xrange(lower, upper))
        
        align_index = np.argmax(scores)
        
        return {
            'shift': align_index + lower,
            'ncc_score': scores[align_index].__float__(),
            'similarity': sim(scores).__float__()
        }
    
    def _align_fast(self, on_data, off_data):
        """
        Internal method to align the <on_data> with the <off_data>
        by shifting in frequencies.
        
        Inputs
        ------
        on_data   : On data matrix.
        off_data  : Off data matrix.
        
        Output
        ------
        The amount the <on_data> needs to shift over the <off_data>.
        
        """
        valid_window = 0 < self.props['align_window'] < on_data.shape[1]
        assert valid_window, 'Error: Invalid Window Size!'
        
        lower = -int(self.props['align_window']/2)
        upper = 1-lower

        X = utils.whiten_data(on_data)
        Y = utils.whiten_data(off_data)
        C = np.dot(X.T, Y)

        scores = np.vectorize(
            lambda i: np.trace(C, offset=i)
        )(range(lower, upper))

        align_index = np.argmax(scores)

        return align_index + lower
    
    def _prmi(self, on_data, off_data, norm=True):
        """
        Internal method to compute the PCA Regional Mutual Information score
        for the given <on_data> and <off_data>.
        
        Inputs
        ------
        on_data   : On data matrix.
        off_data  : Off data matrix.
        
        Output
        ------
        Float value for the PRMI score.
        
        """        
        N, d = on_data.shape
        radius = self.props['prmi_radius']
        
        assert (2 * radius) + 1 < N, 'Error: Selected Radius Too Large!'
        
        Y, X = utils.build_rolling_index(N, d, radius)
        
        on_fv = utils.pca_feature_vector(on_data[Y, X])
        off_fv = utils.pca_feature_vector(off_data[Y, X])
        
        return utils.mutual_info(
            on_fv, 
            off_fv, 
            bins=self.props['histogram_bins'], 
            norm=norm
        )
    
    def _mi(self, on_data, off_data, norm=True):
        """
        Internal method to compute the Mutual Information score for the 
        given <on_data> and <off_data>.
        
        Inputs
        ------
        on_data   : On data matrix.
        off_data  : Off data matrix.
        
        Output
        ------
        Float value for the mutual information score.
        
        """
        return utils.mutual_info(
            on_data.flatten(), 
            off_data.flatten(), 
            bins=self.props['histogram_bins'], 
            norm=norm
        )
    
    def _pcc(self, on_data, off_data):
        """
        Internal method to compute the Pearson Correlation Coefficient scores
        for the given <on_data> and <off_data>. The scores computed are of 4
        variants:
            - global image
            - spatial image
            - distribution of time component
            - distribution of frequency component
        
        Inputs
        ------
        on_data   : On data matrix.
        off_data  : Off data matrix.
        
        Output
        ------
        Dictionary of scores, each is a float value for the PCC score.
        
        """
        on_fv = utils.whiten_data(on_data)
        off_fv = utils.whiten_data(off_data)
        
        on_time_pc, _, on_freq_pc = utils.eigh_extract_component(on_fv, compute_uv=True)
        off_time_pc, _, off_freq_pc = utils.eigh_extract_component(off_fv, compute_uv=True)
        
        on_pc = utils.outer_product(on_time_pc, on_freq_pc)
        off_pc = utils.outer_product(off_time_pc, off_freq_pc)
        
        return {
            'global': utils.pcc_estimator(on_fv, off_fv),
            'time_c': utils.pcc_estimator(on_time_pc, off_time_pc),
            'freq_c': utils.pcc_estimator(on_freq_pc, off_freq_pc),
            'spatial': utils.pcc_estimator(on_pc, off_pc)
        }
    
    def _ssim(self, on_data, off_data):
        """
        Internal method to compute the Structural Similarity Index score
        for the given <on_data> and <off_data>.
        
        Inputs
        ------
        on_data   : On data matrix.
        off_data  : Off data matrix.
        
        Output
        ------
        Float value for the SSIM score.
        
        """
        return compare_ssim(
            utils.remap_data(on_data, max_val=1, dtype='f4'),
            utils.remap_data(off_data, max_val=1, dtype='f4'),
            gaussian_weights=True
        )
        
    def _score_one(self, on_fp, off_fp):
        """
        Internal method to compute the correlator metrics on the given
        .fits pair, <on_fp> and <off_fp> for both the unaligned and
        aligned data.
        
        Inputs
        ------
        on_fp   : On .fits file path.
        off_fp  : Off .fits file path.
        
        Output
        ------
        Dictionary of scores and metrics for the given .fits pair.
        
        """
        on_data = utils.open_fits(on_fp)
        off_data = utils.open_fits(off_fp)
        
        json_obj = { 
            'name': path.basename(on_fp)[:-5], 
            'unaligned': {}, 
            'aligned': {}
        }
        
        for key, score in self.keys:
            json_obj['unaligned'][key] = score(on_data, off_data)
        
        shift = self._align_fast(on_data, off_data)
        json_obj['aligned']['shift'] = shift
        
        if shift != 0:
            on_data, off_data = self._translate(on_data, off_data, shift)
            for key, score in self.keys:
                json_obj['aligned'][key] = score(on_data, off_data)
        else:
            json_obj['aligned'].update(json_obj['unaligned'])
        
        return json_obj
    
    def _score_multi(self, data_dir, show_progress=True):
        """
        Internal method to process the FITS pairs in the given <data_dir>.
        Runs sequentially if multiprocessing is disabled.
        
        Input
        -----
        data_dir  : Path to directory containing .fits files.
        
        Output
        ------
        An array of Python dictionaries, each for every FITS pair.
        
        """
        stdout.write(u'\u250c\u2500 Processing Directory: %s\n' %data_dir)
        stdout.flush()
        bar = utils.progressbar(u'\u2514\u2500 Progress:')
        
        pairs = utils.get_on_off_pairs(data_dir, show_progress=False)
        collection = []
        
        if self.props['multi_process']:
            q = Queue()
            N = len(pairs)
            C = cpu_count()
            chunksize = (N/C) if (N % C == 0) else ((N // C) + 1)

            def _worker(queue, data_slice):
                queue.put([self._score_one(on_fp, off_fp) for on_fp, off_fp in data_slice])

            procs = []
            for i in range(C):
                args = (q, pairs[i * chunksize : (i+1) * chunksize])
                procs.append(Process(target=_worker, args=args))

            for p in procs:
                p.start()

            for _ in bar(range(C)):
                collection.extend(q.get())

            for p in procs:
                p.join()
        else:
            for on_fp, off_fp in bar(pairs):
                collection.append(self._score_one(on_fp, off_fp))
        
        return collection
    
    def _json_writer(self, obj, write_fp):
        """
        Internal method to write Python dictionaries as JSON file.
        
        Inputs
        ------
        obj       : The Python dictionary to jsonify.
        write_fp  : Destination file path.
        
        """
        with open(write_fp, 'w') as fp:
            json.dump(obj, fp, indent=4, separators=(',', ': '))
    
    def run(self):
        """
        Run the correlator and generates the output in the configured output
        directory.
        
        """
        err_msg = 'Error: No Data Pair or Data Directory Assigned!'
        assert self.props['data_pair'] or self.props['data_dir'], err_msg
        
        out_dir = self.props['out_dir']
        
        if not path.exists(out_dir):
            makedirs(out_dir)
        
        if self.props['data_pair']:
            json_obj = self._score_one(*self.props['data_pair'])
            write_fp = path.join(out_dir, json_obj['name'] + '.json')
            self._json_writer(json_obj, write_fp)
        else:
            for i, data_dir in enumerate(utils.get_raw_dirs(self.props['data_dir'])):
                
                json_collection = self._score_multi(data_dir)
                
                write_dir = path.join(out_dir, str(i).zfill(4))
                makedirs(write_dir)
                self._json_writer(
                    {
                        'directory': data_dir,
                        'count': len(json_collection)
                    },
                    path.join(write_dir, '_MANIFEST')
                )
                
                if self.props['coalesce']:
                    write_fp = path.join(write_dir, 'collection.json')
                    self._json_writer(json_collection, write_fp)
                else:
                    for json_obj in json_collection:
                        write_fp = path.join(write_dir, json_obj['name'] + '.json')
                        self._json_writer(json_obj, write_fp)