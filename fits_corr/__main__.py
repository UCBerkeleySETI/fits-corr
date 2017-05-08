from correlator import Correlator
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def _description():
    return '\n'.join(
        [
            u'\u250c\u2500 Correlation & Similarity Analysis for FITS data files.',
            u'\u251c\u2500 Version: 0.1.0-alpha',
            u'\u2514\u2500 \u00a9 Pragaash Ponnusamy 2017'
        ]
    )

def _build_engine(args):
    
    engine = Correlator()
    
    if args.conf:
        engine.load_config(*args.conf)
    else:
        
        if args.dir:
            engine.with_data_dir(*args.dir)
        else:
            engine.with_on_off_pair(*args.file)
            
        bool_operators = [
            (args.g, engine.enable_coalescing),
            (args.m, engine.enable_multiprocessing)
        ]
        
        valued_operators = [
            (args.out, engine.with_output_dir),
            (args.bins, engine.set_histogram_bins),
            (args.prmi, engine.set_prmi_radius),
            (args.w, engine.set_align_window),
            (args.save, engine.save_config)
        ]
        
        for opt, func in bool_operators:
            _ = opt and func()
        
        for opt, func in valued_operators:
            _ = opt and func(opt)
    
    return engine
            

if __name__ == '__main__':

    # Initialize Parser
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=_description()
    )
    
    # Initialize Base Argument Group 
    arg_group = parser.add_mutually_exclusive_group(required=True)
    
    # Base Argument: Config File
    _config = arg_group.add_argument(
        '-c', 
        '--conf', 
        metavar='FILE', 
        help='launch with a .json config file',
        nargs=1
    )
    
    # Base Argument: Input Directory
    _dir = arg_group.add_argument(
        '-d',
        '--dir',
        help='data directory containing .fits files',
        nargs=1
    )
    
    # Base Group: File Pair
    _file = arg_group.add_argument(
        '-f',
        '--file',
        metavar=('ON', 'OFF'),
        nargs=2,
        help='path to on and off file pairs'
    )
    
    # Initialize Options Argument Groups
    opt_group = parser.add_argument_group('options')

    # Options Group: Enable Multiprocessing
    _mp = opt_group.add_argument(
        '-m',
        action='store_true',
        help='enable multiprocessing'
    )
    
    # Options Group: Enable Coalescing
    _coalesce = opt_group.add_argument(
        '-g',
        action='store_true',
        help='enable coalescing'
    )
    
    # Options Group: Output Directory
    _out = opt_group.add_argument(
        '-o',
        '--out',
        metavar='DIR',
        help='path to output directory'
    )
    
    # Options Group: Align Window
    _align = opt_group.add_argument(
        '-w',
        metavar='WIDTH',
        type=int,
        help='alignment window width'
    )
    
    # Options Group: PRMI Radius
    _prmi = opt_group.add_argument(
        '-prmi',
        metavar='RADIUS',
        type=int,
        help='neighborhood window radius for PRMI'
    )
    
    # Options Group: Bins
    _bins = opt_group.add_argument(
        '-bins',
        type=int,
        help='histogram bins for mutual information'
    )
    
    # Options Group: Save
    _save = opt_group.add_argument(
        '-s',
        '--save',
        metavar='PATH',
        help='savepath for config file'
    )

    # Parse Arguments
    engine = _build_engine(parser.parse_args())
    
    # Run Engine
    engine.run()