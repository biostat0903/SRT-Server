import sys
import getopt
import anndata as ad
#
long_opts_list = ['spec=', 'ref=']
param_dict = {'spec': None, 'ref': None}
opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)

for opt, arg in opts:
    if opt == "--spec": 
        param_dict['spec'] = arg
    elif opt == "--ref": 
        param_dict['ref'] = arg
#
spec = param_dict['spec']
ref = param_dict['ref']

exec(open('/public/home/biostat03/project/stwebProject/01_code/srt_server/decon.py').read())
proj_path = "/public/home/biostat03/project/stwebProject/"
ref_path = f'{proj_path}02_data/reference_data/SRT-Server/scRNA-seq/{spec}/{ref}/'

if not os.path.exists(f'{ref_path}cell2location/'):
    os.mkdir(f'{ref_path}cell2location_test/')

adata_sc = ad.read_h5ad(f'{ref_path}/tangram/sc.h5ad')
inf_aver_train(adata_sc=adata_sc, 
               out_path=f'{ref_path}/cell2location/')

