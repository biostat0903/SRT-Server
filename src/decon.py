# deconvolution spot referring to scRNA-seq or marker genes
# up-stream procedure: qc.R
# down-stream procedure: dc_plt.R

import os
import numpy as np
import pandas as pd
import anndata as ad
from h5py import File
from scipy.sparse import csr_matrix
from re import sub

# method_path setting
refer_path = '/public/home/biostat03/project/stwebProject/02_data/reference_data/SRT-Server/scRNA-seq/'

# Fix parematers
NUM_EPOCH_TG = 1000
NUM_EPOCH_C2L_TRAINING = 10
NUM_EPOCH_C2L_FIT = 100
NUM_EPOCH_C2L_FIT = 10
SC_MAP_TRESH = 0.01                ## thresholod of elements in the sc-mapping matrix

# Function 1: train single cell reference model for cell2location
def inf_aver_train(adata_sc, 
                   out_path: str=None
):
    import cell2location
    ## filter by genes
    adata_sc.var['gene'] = adata_sc.var.index
    selected = cell2location.utils.filtering.filter_genes(adata_sc,
                                                          cell_count_cutoff=5,
                                                          cell_percentage_cutoff2=0.03,
                                                          nonz_mean_cutoff=1.12)
    ## filter the object
    adata_sc = adata_sc[:, selected].copy()
    ## fit regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_sc,
                                                       batch_key='sampleInfo',
                                                       labels_key='cellType')
    ## create the regression model
    mod = cell2location.models.RegressionModel(adata_sc)
    mod.train(max_epochs=NUM_EPOCH_C2L_TRAINING,
              batch_size=2500,
              train_size=1,
              lr=0.002)
    ## save model and sc-anndata object with results
    mod.save(out_path, overwrite=True)
    ## export the estimated cell abundance (summary of the posterior distribution).
    adata_sc = mod.export_posterior(adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    ## export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_sc.uns['mod']['factor_names']
    ## output
    out_file = f'{out_path}/inf_aver.txt.gz'
    inf_aver.to_csv(out_file, header=True, index=True,
                       sep='\t', mode='w', compression='gzip')
    return out_file

# Fucntion 2: transformat to single cell data
def sc_format(inpath, 
              outpath,
              decon_method: str=None
):
    if (os.path.exists(f'{inpath}/sc_count.txt') and os.path.exists(f'{inpath}/sc_meta.txt')):
        sc_count = pd.read_csv(f'{inpath}/sc_count.txt', sep = "\t")
        sc_meta = pd.read_csv(f'{inpath}/sc_meta.txt', sep = "\t")
    else:
        raise SystemExit('Data format of single cell reference is invalid!')
    intersect_cell = np.intersect1d(sc_count.columns, sc_meta.index)
    sc_count_inter = sc_count.loc[:, intersect_cell].copy()
    sc_meta_inter = sc_meta.loc[intersect_cell, ["cellType","sampleInfo"]].copy()
    sc_ref = ad.AnnData(X = sc_count_inter.T,
                        obs=pd.DataFrame(sc_meta_inter, index=intersect_cell),
                        var=pd.DataFrame(sc_count_inter.index, columns=['gene'], index=sc_count_inter.index))
    sc_ref.var_names.name = "gene"
    sc_ref.obs_names.name = "cell"
    if decon_method == 'cell2location':
        sc_ref_t = inf_aver_train(sc_ref, outpath)
        return outpath
    else:
        ref_adata_file = f"{outpath}/sc.h5ad"
        sc_ref.write(ref_adata_file)
        return outpath

# Function 3: check function
def decon_check(data_path: str = None,                ## String: output path for qc procedure
                source_panel: str = None,             ## String: the source of panel: {Local} {Cloud}
                ref_path: str = None,                 ## String: reference panel path
                decon_method = 'cell2location',       ## String: deconvolution methods: {cell2location} and {Tangram}
                species: str = None,                  ## String: species: {Human} and {Mouse}
                status: str = None,                   ## String: combine tissue and disease status
                c2l_cells_per_loc = 30,               ## Int: cell2loaction: cell per location
                c2l_det_alpha = 200,                  ## Int: cell2loaction: detection alpha
                tg_mode = 'clusters',                 ## String: Tangram
                scMapping = False,                    ## Boolean: deconvolution methods: {cell2location} and {Tangram}
                out_path: str = None                  ## String: output path for decon procedure
):
    ## check qc file
    check_file = f'{data_path}/qc_call_file.txt'
    if not os.path.exists(check_file):
        raise SystemExit('No QC file! Please run QC module!')
    else:
        qc_param = pd.read_csv(check_file, header=None)
        spatial_data_filename = qc_param.at[0,0]
        sample_size = qc_param.at[1,0]
    ## check reference
    if source_panel == 'Local':
        ref_out_dir = f'{out_path}/ref_out/'
        if not os.path.exists(ref_out_dir):
            os.makedirs(ref_out_dir, exist_ok=True)
        sc_dir = sc_format(ref_path, ref_out_dir, decon_method)
    else:
        sc_dir = f'{refer_path}/{species}/{status}/{decon_method}/'
    ## output 
    check_file_path = f'{out_path}/decon_check_file.txt'
    check_file = pd.DataFrame([spatial_data_filename, sample_size, sc_dir, decon_method, 
                               c2l_cells_per_loc, c2l_det_alpha, tg_mode, scMapping, source_panel])
    check_file.to_csv(check_file_path,
                      header=None, index=None, sep='\t', mode='w')
    return 0

# Function 4: create anndata object using h5 file
def h5_2_adata(h5_dat,
               sample_idx, 
               image=False,
):
    # extract count
    count_s = h5_dat['count'][f'count.s{sample_idx}'][:]
    cellID_s = h5_dat['cellID'][f'cellID.s{sample_idx}'][:]
    featureID_s = h5_dat['featureID'][f'featureID.s{sample_idx}']
    # remove bytes code
    cellID_s = [bytes.decode(idx) for idx in cellID_s]
    featureID_s = [bytes.decode(idx) for idx in featureID_s]
    # extract coord
    coord_s = h5_dat['meta.data'][f'meta.data.s{sample_idx}'][:].T
    # make adata
    adata_stx = ad.AnnData(X=count_s,
                          obs=pd.DataFrame(abs(coord_s), columns=['x', 'y'], index=cellID_s),
                          obsm={"spatial": abs(coord_s),
                                "orig_coord": pd.DataFrame(coord_s, columns=['x', 'y'], index=cellID_s)},
                          var=pd.DataFrame(featureID_s, columns=['gene'], index=featureID_s))
    adata_stx.var_names.name = "gene"
    adata_stx.obs_names.name = "cell"
    adata_stx.obs['sample'] = sample_idx
    if image == True:
        hires_img = h5_dat['image.data'][f'image.data.s{sample_idx}'][:]
        image_meta = h5_dat['image.meta'][f'image.meta.s{sample_idx}'][:]
        image_meta = [bytes.decode(idx) for idx in image_meta]
        scalefactors = {"fiducial_diameter_fullres": float(image_meta[1]),
                        "spot_diameter_fullres": float(image_meta[2]),
                        "tissue_hires_scalef": float(image_meta[3]),
                        "tissue_lowres_scalef": float(image_meta[4])}
        adata_stx.uns["spatial"] = {f's{sample_idx}': {"image": {"hires": hires_img},
                                                       "image_path": image_meta[0],
                                                       "scalefactors": scalefactors}}
    return adata_stx

# Function 5: convert abundance/prob to proprotion
def convert_density_to_prop(df):
    rowSum = df.sum(axis = 1)
    prop = df.div(rowSum, axis = 0)
    return prop

# Function 6: cell2location call function
def cell2loc_call(h5_file: str=None,
                  ref_path: str=None,
                  out_path: str=None,
                  sc_mapping=False,
                  sample_size=1,
                  cells_per_loc=30,
                  det_alpha=200
):
    ## import extra modules
    import cell2location
    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes
    
    ## load data
    print("Loading data...")
    h5_dat = File(h5_file, "r")
    ### format st adata
    if sample_size == 1:
        adata_st = h5_2_adata(h5_dat, sample_size)
    else:
        ### create each anndata object from h5 file
        slides = []
        for sample_idx in range(1,sample_size+1):
            slides.append(h5_2_adata(h5_dat, sample_idx))
        ### combine anndata objects together
        adata_st = slides[0].concatenate(
            slides[1:],
            batch_key="idx",
            uns_merge="unique",
            index_unique=None)
    ## load training result
    inf_aver = pd.read_csv(f'{ref_path}/inf_aver.txt.gz', 
                           sep = "\t", header = 0, index_col = 0)
    ## find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_st.var_names, inf_aver.index)
    adata_st = adata_st[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    ## fit cell2location
    ### prepare anndata for cell2location model
    print("Fitting cell2location on spatial data...")
    cell2location.models.Cell2location.setup_anndata(adata=adata_st, batch_key="sample")
    map_path = f"{out_path}/cell2location_map"
    ### create and train the model
    if os.path.exists(f"{map_path}/model.pt"):
        mod = cell2location.models.Cell2location.load(f"{map_path}/", adata_st)
    else:
        ### N_cells_per_location: the expected average cell abundance: tissue-dependent ######
        ###                       hyper-prior which can be estimated from paired histology. ##
        ### detection_alpha: hyperparameter controlling normalisation of #####################
        ###                  within-experiment variation in RNA detection ####################
        #
        mod = cell2location.models.Cell2location(adata_st,
                                                 cell_state_df=inf_aver,
                                                 N_cells_per_location=cells_per_loc,
                                                 detection_alpha=det_alpha)
        #
        mod.train(max_epochs=NUM_EPOCH_C2L_FIT,
                  batch_size=None,
                  train_size=1)
        ### Save model and st-anndata object with results
        map_path = f"{out_path}/cell2location_map/"
        if not os.path.exists(map_path):
            os.mkdir(map_path)
        
        mod.save(map_path, overwrite=True)
    
    ### Export the estimated cell abundance (summary of the posterior distribution).
    adata_st = mod.export_posterior(
        adata_st, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )
    
    ## process output 
    if sc_mapping == str(True):
        # map file (sc mapping)
        ## Compute expected expression per cell type
        adata_st.obs[adata_st.uns['mod']['factor_names']] = adata_st.obsm['q05_cell_abundance_w_sf']
        expected_dict = mod.module.model.compute_expected_per_cell_type(
            mod.samples["post_sample_q05"], mod.adata_manager
        )
        ## Add to anndata layers
        for i, n in enumerate(mod.factor_names_):
            adata_st.layers[n] = expected_dict['mu'][i]
        #
        all_ct = list(adata_st.layers.keys())
        for ctx in all_ct:
            adata_st.layers[ctx][adata_st.layers[ctx] < SC_MAP_TRESH] = 0
    
    ## format deconv output
    adata_st.obs[adata_st.uns['mod']['factor_names']] = convert_density_to_prop(adata_st.obsm['q05_cell_abundance_w_sf'])
    adata_st.obs['cell'] = adata_st.obs_names
    cols = adata_st.obs.columns.tolist()
    adata_st.obs = adata_st.obs.loc[:,cols[2:3]+cols[-1:]+cols[6:-1]]
    
    # output data frame
    print("Saving results...")
    st_adata_file = f"{out_path}/c2l_map_sp.h5ad"
    adata_st.write(st_adata_file)
    
    return st_adata_file

# Function 7: image process for tangram_constrained 
def tg_img_process(adata_st, scalef, 
                   img_file: str=None
):
    # import extra module
    from squidpy import im
    # process image
    img = im.ImageContainer(img=img_file,
                            layer='image',
                            lazy=True,
                            scale=scalef)
    im.process(img=img,
               layer="image",
               method="smooth")
    im.segment(img=img,
               layer="image_smooth",
               method="watershed",
               channel=0)
    # define image layer to use for segmentation
    features_kwargs = {
        "segmentation": {
            "label_layer": "segmented_watershed",
            "props": ["label", "centroid"],
            "channels": [1, 2],
        }
    }
    # calculate segmentation features
    im.calculate_image_features(adata_st,
                                img,
                                layer="image",
                                key_added="image_features",
                                features_kwargs=features_kwargs,
                                features="segmentation",
                                mask_circle=True)
    
    # determine cell count
    adata_st.obs["cell_count"] = adata_st.obsm["image_features"]["segmentation_label"]
    
    return adata_st

# Function 8: tangram call function
def tg_call(h5_file: str=None,
            ref_path: str=None,
            tg_mode="clusters",
            out_path: str=None,
            sample_idx=1,
            sc_mapping=False
):
    ## import extra modules
    import tangram as tg
    ## load data
    print("Loading data...")
    adata_sc = ad.read_h5ad(f'{ref_path}/sc.h5ad')
    h5_dat = File(h5_file, "r")
    if tg_mode == "constrained":
        image=True
    else:
        image=False
    
    adata_st = h5_2_adata(h5_dat, sample_idx, image=image)
    
    ## Estimate how many cells are present in each voxel
    if tg_mode == "constrained":
        print("Loading image information for constrained mode……")
        scalef_hires = adata_st.uns['spatial'][f's{sample_idx}']['scalefactors']['tissue_hires_scalef']
        img_file = adata_st.uns['spatial'][f's{sample_idx}']['image_path'][:]
        adata_st = tg_img_process(adata_st=adata_st,
                                  img_file=img_file,
                                  scalef=scalef_hires)
        target_count = adata_st.obs.cell_count.sum()
        density_prior = np.array(adata_st.obs.cell_count) / adata_st.obs.cell_count.sum()
    else:
        target_count = None
        density_prior = "rna_count_based"
    
    ## Deconvolution via alignment
    print("Run tangram mapping……")
    tg.pp_adatas(adata_sc, adata_st)
    #
    ad_map = tg.map_cells_to_space(adata_sc,
                                   adata_st,
                                   mode=tg_mode,
                                   cluster_label="cellType",
                                   target_count=target_count,
                                   density_prior=density_prior,
                                   num_epochs=NUM_EPOCH_TG,
                                   device='cpu')
    
    ## Export the mapping between spot ID and segmentation ID to adata.uns.
    if tg_mode == "constrained":
        tg.create_segment_cell_df(adata_st)
    
    ## Map cell types as result of the deconvolution step to putative segmentation ID.
    tg.project_cell_annotations(ad_map,
                                adata_st,
                                annotation="cellType")
    # format deconv output
    adata_st.obs = convert_density_to_prop(adata_st.obsm['tangram_ct_pred'])
    adata_st.obs['sample'] = sample_idx
    adata_st.obs['cell'] = adata_st.obs_names
    cols = adata_st.obs.columns.tolist()
    adata_st.obs = adata_st.obs.loc[:,cols[-2:]+cols[:-2]]
    ## map file (sc mapping)
    if sc_mapping == str(True) and tg_mode != "clusters":
        # calculate cell type-specific expression for each spot
        all_ct = np.unique(ad_map.obs['cellType'])
        all_spot = ad_map.var_names
        gene_op = adata_sc.uns['overlap_genes']
        adata_st.obs = adata_st[:, gene_op].obs.copy()
        adata_st.obsm = {"orig_coord": adata_st.obsm['orig_coord']}
        gene_op_idx = np.where(np.isin(adata_sc.var_names, gene_op))[0]
        sc_mat_map = csr_matrix(adata_sc.X[:,gene_op_idx])
        for ctx in all_ct[12:]:
            prob_spot = ad_map.X.copy()
            for spotx in all_spot:
                spot_idx = np.where(ad_map.var_names == spotx)[0]
                # prob_spot = ad_map.X[:,spot_idx]
                ct_idx = np.where(ad_map.obs['cellType'] != ctx)[0]
                prob_spot[ct_idx,spot_idx] = 0
            
            adata_st.layers[ctx] = csr_matrix(prob_spot.T@sc_mat_map)
    ## output data frame
    print("Saving results...")
    st_adata_file = f"{out_path}/tg_{tg_mode}_map_sp_s{sample_idx}.h5ad"
    adata_st.write(st_adata_file)
    
    return st_adata_file

# Function 9: format and output sc mapping results
def out_mapping(adata_st_sample, 
                out_path_sample: str=None
):
    # combine sc mapping results
    all_ct = list(adata_st_sample.layers.keys())
    scmap_count = pd.concat([pd.DataFrame(adata_st_sample.layers[ctx].toarray(),
                                          columns = adata_st_sample.var_names) for ctx in all_ct])
    scmap_coord = pd.concat([adata_st_sample.obsm['orig_coord'] for ctx in all_ct])
    scmap_meta = pd.concat([pd.DataFrame({"cell": [spotx+"_"+ctx for spotx in adata_st_sample.obs_names],
                                          "celltype": ctx}) for ctx in all_ct])
    # write out
    scmap_count_file = f'{out_path_sample}/count.txt.gz'
    scmap_coord_file = f'{out_path_sample}/coord.txt.gz'
    scmap_meta_file = f'{out_path_sample}/annotation.txt.gz'
    scmap_count.to_csv(scmap_count_file, header=True, index=None,
                       sep='\t', mode='w', compression='gzip')
    scmap_coord.to_csv(scmap_coord_file, header=True, index=None,
                       sep='\t', mode='w', compression='gzip')
    scmap_meta.to_csv(scmap_meta_file, header=True, index=None,
                      sep='\t', mode='w', compression='gzip')
    
    return str(scmap_count_file+","+scmap_coord_file+","+scmap_meta_file)

# Function 10: call function
def decon_call(out_path: str=None,                     ## String: output path for decon procedure
):
    ## set parameters
    decon_check = pd.read_csv(f'{out_path}/decon_check_file.txt', header=None)
    spatial_data_filename = decon_check.at[0,0]
    sample_size = int(decon_check.at[1,0])
    sc_dir = decon_check.at[2,0]
    decon_method = decon_check.at[3,0]
    c2l_cells_per_loc = int(decon_check.at[4,0])
    c2l_det_alpha = int(decon_check.at[5,0])
    tg_mode = decon_check.at[6,0]
    sc_mapping = decon_check.at[7,0]
    source_panel = decon_check.at[8,0]
    ## two methods
    if decon_method == 'tangram':
        print(f"Run tangram in {tg_mode} mode……")
        call_file = []
        for sample_idx in list(range(1, sample_size+1)):
            call_filex = tg_call(h5_file = spatial_data_filename,
                                 ref_path = sc_dir,
                                 sc_mapping = sc_mapping,
                                 sample_idx = sample_idx,
                                 tg_mode = tg_mode,
                                 out_path = out_path)
            call_file.append(call_filex)
    else:
        print("Run cell2location……")
        call_filex = cell2loc_call(h5_file = spatial_data_filename,
                                   ref_path = sc_dir,
                                   sc_mapping = sc_mapping,
                                   sample_size = sample_size,
                                   cells_per_loc = c2l_cells_per_loc,
                                   det_alpha = c2l_det_alpha,
                                   out_path = out_path)
        call_file = [call_filex]
    ## output call file
    print("Saving call file……")
    call_file_path = f"{out_path}/decon_call_file.txt"
    call_file = pd.DataFrame(call_file)
    call_file.to_csv(call_file_path,
                     header=None, index=None, sep='\t', mode='w')
    ## delete temportary files
    if source_panel == 'Local':
        os.system(f'rm -rf {sc_dir}')
    return 0

# Function 11: post function
def decon_post(out_path: str=None                     ## String: output path for decon procedure
):
    ## load data
    print("Post processing……")
    decon_check = pd.read_csv(f'{out_path}/decon_check_file.txt', header=None)
    decon_call = pd.read_csv(f'{out_path}/decon_call_file.txt', header=None)
    sample_size = int(decon_check.at[1,0])
    decon_method = decon_check.at[3,0]
    sc_mapping = decon_check.at[7,0]
    result_dir = f'{out_path}/decon_result/{decon_method}/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir, exist_ok=True)
    ## set parameter for different methods
    if decon_method == "cell2location":
        adata_st = ad.read_h5ad(decon_call.iloc[0,0])
        dc_method_mode = decon_method
    if decon_method == "tangram":
        tg_mode = decon_check.at[6,0]
        dc_method_mode = f'{decon_method}_{tg_mode}'
    ## format output
    decon_file = f'{result_dir}/decon_Matrix.txt'
    post_file = [decon_method, decon_file]
    decon_list = []
    for sample_idx in list(range(1, sample_size+1)):
        ### set adata_st for sample_idx
        if decon_method == "tangram":
            adata_stx = ad.read_h5ad(decon_call.iloc[sample_idx-1,0])
        else:
            sample_cell_idx = np.where(adata_st.obs['sample'] == sample_idx)[0]
            adata_stx = adata_st[sample_cell_idx, :].copy()
        ### main results
        adata_stx.obs.columns = [sub(r'\s', "_", namex) for namex in adata_stx.obs.columns]
        decon_list.append(adata_stx.obs)
        ### sc mapping results
        if sc_mapping == str(True):
            out_path_samplex = f'{result_dir}/sample{sample_idx}/'
            if not os.path.exists(out_path_samplex):
                os.mkdir(out_path_samplex)
            scmap_pathx = out_mapping(adata_st_sample=adata_stx,
                                      out_path_sample=out_path_samplex)
            post_file.append(scmap_pathx)
    decon_df = pd.concat(decon_list)
    decon_df.to_csv(decon_file, header=True, index=None, sep='\t', mode='w')
    ## post file
    decon_post_path = pd.DataFrame(post_file)
    decon_post_path.to_csv(f'{out_path}/decon_post_file.txt',
                           header=None, index=None, sep='\t', mode='w')
    return 0
