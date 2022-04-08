import argparse
import numpy as np
from multiprocessing import Queue,Process
from tqdm import tqdm

from quality_net_utilities import *
from segment_slide_wbc_rbc import *
from characterize_cells import *
from cells_to_annotations import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Predicts which tiles are of good quality in WBS.')

    # slide arguments
    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',type=str,default=None,
                        help="Path to slide.")
    parser.add_argument('--input_height',dest = 'input_height',
                        action = 'store',type = int,default = 512,
                        help = 'The file extension for all images.')
    parser.add_argument('--input_width',dest = 'input_width',
                        action = 'store',type = int,default = 512,
                        help = 'The file extension for all images.')
    parser.add_argument('--rescale_factor',dest='rescale_factor',
                        action='store',type=float,default=1,
                        help='Factor to resize cells (due to different mpp).')

    # qc arguments
    parser.add_argument('--checkpoint_path_qc',dest = 'checkpoint_path_qc',
                        action = 'store',type = str,default = 'summaries',
                        help = 'Path to checkpoint.')
    parser.add_argument('--batch_size_qc',dest = 'batch_size_qc',
                        action = 'store',type = int,default = 4,
                        help = 'Size of mini batch for QC.')
    parser.add_argument('--output_qc',dest = 'output_qc',
                        action = 'store',type = int,default = 4,
                        help = 'Output for QC.')

    # segmentation arguments
    parser.add_argument('--unet_checkpoint_path',dest='unet_checkpoint_path',
                        action='store',default=None,
                        help='Path to U-Net checkpoint.')
    parser.add_argument('--depth_mult',dest='depth_mult',
                        action='store',type=float,default=1.0,
                        help='Depth of the U-Net.')
    parser.add_argument('--wbc_output_path_seg',dest='wbc_output_path_seg',
                        action='store',default=None,
                        help='Output path for the WBC segmentation file.')
    parser.add_argument('--rbc_output_path_seg',dest='rbc_output_path_seg',
                        action='store',default=None,
                        help='Output path for the RBC segmentation file.')

    # characterization arguments
    parser.add_argument('--output_path_wbc_char',dest='output_path_wbc_char',
                        default=None,action='store',
                        help='Output path for the WBC characterization file.')
    parser.add_argument('--output_path_rbc_char',dest='output_path_wbc_char',
                        default=None,action='store',
                        help='Output path for the WBC characterization file.')
    parser.add_argument('--n_processes',dest='n_processes',
                        action='store',default=2,type=int,
                        help='Number of concurrent processes for characterization')
    parser.add_argument('--standardizer_params',dest='standardizer_params',
                        action='store',default="scripts/python/rbc_scaler_params",
                        type=str,help='Path to standardizer parameters (RBC char).')
    parser.add_argument('--xgboost_model',dest='xgboost_model',
                        action='store',default="scripts/python/rbc_xgb_model",
                        type=str,help='Path to XGBoost model (RBC char only).')

    args = parser.parse_args()

    # define step-specific networks and processes
    ## define QC network
    quality_net = keras.models.load_model(args.checkpoint_path)
    def generator():
        G = image_generator_slide(
            args.slide_path,args.input_height,args.input_width)
        for image,coords in G:
            image = image / 255.
            yield image,coords
    output_types = (tf.float32,tf.string)
    output_shapes = (
        [args.input_height,args.input_width,n_channels],[])
    tf_dataset = tf.data.Dataset.from_generator(
        generator,output_types=output_types,output_shapes=output_shapes)
    tf_dataset = tf_dataset.batch(
        args.batch_size_qc,
        drop_remainder=False)
    tf_dataset = tf_dataset.prefetch(5)

    ## define segmentation network and processes
    h,w,extra = args.input_height,args.input_width,128
    new_h,new_w = h+extra,w+extra
    u_net = UNet(depth_mult=args.depth_mult,padding='SAME',
                factorization=False,n_classes=2,
                dropout_rate=0,squeeze_and_excite=False)
    u_net.load_weights(args.unet_checkpoint_path)
    u_net.trainable = False
    u_net.make_predict_function()

    wbc_process = WBCProcess(args.wbc_output_path_seg,h,w,extra)
    wbc_process.init_hdf5()
    rbc_process = RBCProcess(args.rbc_output_path_seg,h,w,extra)
    rbc_process.init_hdf5()

    ## characterization processes
    feature_extractor_wbc = FeatureExtractor(args.rescale_factor,"wbc")
    wbc_queue = Queue()
    rbc_results_queue = Queue()
    def wbc_char_fn(qqo):
        q,qo = qqo
        proceed = True
        while proceed == True:
            x = wbc_queue.get()
            if x != None:
                x = image,mask,cnt,C
                f = feature_extractor_rbc.extract_features(
                    image,mask,cnt)
                qo.put([f,C])
            else: proceed = False
    qo.put(None)
    wbc_feature_process = Process(
        wbc_char_fn,[wbc_queue,wbc_results_queue])

    feature_extractor_rbc = FeatureExtractor(args.rescale_factor,"rbc")
    rbc_queue = Queue()
    rbc_results_queue = Queue()
    def rbc_char_fn(qqo):
        q,qo = qqo
        proceed = True
        while proceed == True:
            x = rbc_queue.get()
            if x != None:
                x = image,mask,cnt,C
                f = feature_extractor_wbc.extract_features(
                    image,mask,cnt)
                qo.put([f,C])
            else: proceed = False
        qo.put(None)
    rbc_feature_process = Process(
        rbc_char_fn,[q,rbc_queue,rbc_results_queue])

    # create output objects
    qc_out = open(args.output_qc,'w')

    for images,coords in tqdm(tf_dataset):
        prediction = quality_net(images)
        for c,p,image in zip(coords.numpy(),prediction.numpy(),image.numpy()):
            out_str = 'OUT,{},{},{}\n'.format(c.decode(),int(p>0.5),float(p))
            qc_out.write(out_str)
            if p > 0.5:
                c = list(c)
                masked_wbc = mask_wbc(n_i,u_net,tta=True)
                image = image.numpy()
                image = np.uint8(image*255)
                
                wbc_process.process_element(
                    [image,masked_wbc,args.rescale_factor,c])
                rbc_process.process_element(
                    [image,args.rescale_factor,c])

    # stopping processes  
    rbc_queue.put(None)
    wbc_queue.put(None)

    # aggregating features for RBC
    all_rbc_features = []
    rbc_cell_centers = []
    all_ks_rbc = []
    getting_rbc = True
    cell_idx = 0
    while getting_rbc == True:
        f = rbc_results_queue.get()
        if f is not None:
            f,C = f
            if np.any(np.isnan(f)) or np.any(np.isinf(f)):
                pass
            else:
                all_ks_rbc.append(cell_idx)
                all_rbc_features.append(f)
                rbc_cell_centers.append(C)
        cell_idx = 0
        else:
            getting_rbc = False
    all_rbc_features = np.array(all_rbc_features)

    standardized_scaler = StandardScaler()
    fitted_model = xgboost.XGBClassifier(
        booster='dart',eval_metric="logloss")
    with open(args.standardizer_params,"rb") as o:
        sc_params = pickle.load(o)

    standardized_scaler.n_features_in_ = len(sc_params['mean'])
    standardized_scaler.mean_ = sc_params['mean']
    standardized_scaler.scale_ = sc_params['std']
    standardized_scaler.var_ = sc_params['std']**2
    fitted_model.load_model(args.xgboost_model)

    if all_rbc_features.size > 0:
        all_rbc_features_ = standardized_scaler.transform(all_rbc_features)
        predictions = fitted_model.predict(all_rbc_features_)
        all_rbc_features = all_rbc_features[predictions,:]
        all_ks = all_ks[predictions]

    with h5py.File(args.output_path_rbc_char,'w') as F_out:
        F_out['cells/0'] = all_rbc_features
        F_out['means'] = np.mean(all_rbc_features,axis=0)
        F_out['variances'] = np.var(all_rbc_features,axis=0)
        F_out['cell_centers'] = np.array(cell_centers,dtype="S")
        F_out['cell_center_idxs'] = all_ks

    # aggregating features for WBC
    all_wbc_features = []
    wbc_cell_centers = []
    all_ks_wbc = []
    getting_wbc = True
    cell_idx = 0
    while getting_wbc == True:
        f = wbc_results_queue.get()
        if f is not None:
            f,C = f
            if np.any(np.isnan(f)) or np.any(np.isinf(f)):
                pass
            else:
                all_ks_wbc.append(cell_idx)
                all_wbc_features.append(f)
                wbc_cell_centers.append(C)
        cell_idx = 0
        else:
            getting_wbc = False
    all_wbc_features = np.array(all_wbc_features)

    with h5py.File(args.output_path_wbc_char,'w') as F_out:
        F_out['cells/0'] = all_rbc_features
        F_out['means'] = np.mean(all_rbc_features,axis=0)
        F_out['variances'] = np.var(all_rbc_features,axis=0)
        F_out['cell_centers'] = np.array(cell_centers,dtype="S")
        F_out['cell_center_idxs'] = all_ks
