from glob import glob
import os
import h5py
import numpy as np
from scipy.spatial import distance
from math import inf

from skimage.transform import rotate
from skimage import io
from PIL import Image
import cv2

def read_image(path):
    return io.imread(path)

def rotate_coords(x,y,angle):
    angle = np.radians(angle)
    new_x = x*np.cos(angle) - y*np.sin(angle)
    new_y = y*np.cos(angle) + x*np.sin(angle)
    return new_x,new_y

def show_image_mask(image,mask=None):
    if mask is not None:
        image = image * np.where(mask > 0,1,0.3)
        image = image.astype(np.uint8)
    return image

def show_image_boxes(image,bounding_boxes):
    for box in bounding_boxes:
        x1,y1,x2,y2 = np.round(box).astype(np.int32)
        cv2.rectangle(image, (y1, x1), (y2, x2), (255,0,0), 1)
    return image

def cart2pol(x,y):
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    return rho,phi

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y

def show_image_polygons(image,bounding_polygons):
    for polygon in bounding_polygons:
        cv2.polylines(image,[polygon.T[:,np.newaxis,:]],True,(255,0,0),
                      thickness=1)
        cv2.polylines(image,polygon.T[:,np.newaxis,:],True,(0,0,255),
                      thickness=1)
    return image

def show_image_polygons_centers(image,bounding_polygons,centers):
    for polygon,center in zip(bounding_polygons,centers):
        cv2.polylines(image,[polygon.T[:,np.newaxis,:]],True,(255,0,0),
                      thickness=1)
        cv2.polylines(image,polygon.T[:,np.newaxis,:],True,(0,0,255),
                      thickness=1)
        for point in polygon.T:

            cv2.line(image,
                     (point[0],point[1],),
                     (center[0],center[1]),(0,255,0))
    return image

def box_prediction_format(size,anchor,obj,scale):
    """
    Formats box predictions from YOLOv3.
    """
    _,x,y = np.where(obj > 0)
    if len(x) > 0:
        centers = size[:2,x,y]
        c_x,c_y = centers[0,:]+x,centers[1,:]+y
        c_x,c_y = (
            np.round(scale[0]*c_x),
            np.round(scale[1]*c_y))
        a_x,a_y = anchor[0,x,y],anchor[1,x,y]
        h,w = (
            np.round(np.exp(size[2,x,y])*a_x*scale[0]),
            np.round(np.exp(size[3,x,y])*a_y*scale[1]))
        bboxes = []
        for a,b,c,d in zip(c_x,c_y,h,w):
            bboxes.append(np.array([
                a - c/2,
                b - d/2,
                a + c/2,
                b + d/2
            ]))
        return bboxes
    else:
        return []

def select_region(record,dimensions):

    x1,y1,x2,y2 = [int(x) for x in dimensions]

    if record['centers'].shape[0] > 0:
        conditions = np.array([
            record['centers'][:,0] > y1,
            record['centers'][:,1] > x1,
            record['centers'][:,0] < y2-1,
            record['centers'][:,1] < x2-1
        ])
        center_idxs = np.where(np.all(conditions,axis=0))[0]
    else:
        center_idxs = np.array([])

    if center_idxs.shape[0] > 0:
        centers = np.int32(record['centers'][center_idxs,:])
        centers[:,0] -= y1
        centers[:,1] -= x1
        bboxes = np.int32(record['bounding_boxes'][center_idxs,:])
        bboxes[:,0] -= x1
        bboxes[:,2] -= x1
        bboxes[:,1] -= y1
        bboxes[:,3] -= y1
        bounding_polygons = {}
        for key in record['bounding_polygons']:
            bp = np.int32(record['bounding_polygons'][key][center_idxs,:])
            bp[:,0,:] -= y1
            bp[:,1,:] -= x1
            bounding_polygons[key] = bp
        edges = {}
        for i,x in enumerate(center_idxs):
            edge = record['edges'][str(x+1)][:,:]
            edge[0,] -= y1
            edge[1,] -= x1
            edges[str(i+1)] = edge
    else:
        centers = np.array([])
        bboxes = np.array([])
        bounding_polygons = {
            x:np.array([]) for x in record['bounding_polygons']}
        edges = {}

    return {
        'image':record['image'][x1:x2,y1:y2],
        'mask': record['mask'][x1:x2,y1:y2],
        'weight_map':record['weight_map'][x1:x2,y1:y2],
        'edges':edges,
        'centers':centers,
        'bounding_boxes':bboxes,
        'bounding_polygons':bounding_polygons
    }

def rotate_image(image,angle,interpolation=cv2.INTER_LINEAR):
    image_center = tuple(np.array(image.shape[1::-1]) / 2)
    rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv2.warpAffine(image,rot_mat,
                            image.shape[1::-1],
                            flags=interpolation)
    return result

def record_rotation(record,angle):
    sh = record['image'].shape
    angle = angle % 360
    record = {
        'image':record['image'][:],
        'mask':record['mask'][:],
        'weight_map':record['weight_map'][:],
        'bounding_boxes':record['bounding_boxes'][:],
        'bounding_polygons':{x:record['bounding_polygons'][x][:]
                             for x in record['bounding_polygons']},
        'centers':record['centers'][:],
        'edges':{x:record['edges'][x][:]
                 for x in record['edges']}
    }

    record['image'] = rotate_image(
        record['image'],angle)
    record['mask'] = rotate_image(
        record['mask'],angle,cv2.INTER_NEAREST)[:,:,np.newaxis]
    record['weight_map'] = rotate_image(
        record['weight_map'],angle)[:,:,np.newaxis]
    if record['centers'].shape[0] > 0:
        record['centers'][:,0],record['centers'][:,1] = rotate_coords(
            record['centers'][:,0]-sh[0]//2,
            record['centers'][:,1]-sh[1]//2,
            -angle
        )
        record['centers'][:,0] += sh[0]//2
        record['centers'][:,1] += sh[1]//2
        centers_out = []
        edges = record['edges']
        edges_out = {}

        bboxes = []
        bounding_polygons = {
            x:[] for x in record['bounding_polygons']
        }
        i = 0
        for edge_key,center in zip(edges,record['centers']):
            if (1<=center[0]<(sh[0]-1)) and (1<=center[1]<(sh[1]-1)):
                i += 1
                centers_out.append(center)
                edge = edges[edge_key]
                edge_x,edge_y = edge[0,:],edge[1,:]
                edge_x,edge_y = rotate_coords(
                    edge_x-sh[0]/2,edge_y-sh[1]/2,-angle)
                edge_x,edge_y = edge_x + sh[0]/2,edge_y + sh[1]/2

                bboxes.append([edge_y.min(),edge_x.min(),
                               edge_y.max(),edge_x.max()])

                edges_out[str(i)] = np.array([edge_x,edge_y])

                for key in bounding_polygons:
                    n = int(key)
                    bounding_polygon = edges_to_polygon(
                        edge_x-center[0],edge_y-center[1],n)

                    bounding_polygon[0] = np.roll(
                        bounding_polygon[0],3*(n//4)-1)
                    bounding_polygon[1] = np.roll(
                        bounding_polygon[1],3*(n//4)-1)

                    bounding_polygon[0] = np.flip(bounding_polygon[0])
                    bounding_polygon[1] = np.flip(bounding_polygon[1])

                    bounding_polygon = np.array(
                        [bounding_polygon[1]+center[0],
                         bounding_polygon[0]+center[1]]
                    )


                    bounding_polygons[key].append(bounding_polygon)

        record['edges'] = edges_out
        record['bounding_boxes'] = np.array(bboxes)
        record['bounding_polygons'] = {
            x:np.array(bounding_polygons[x])
            for x in bounding_polygons}
        record['centers'] = np.array(centers_out)

    return record

def edges_to_polygon(edge_x,edge_y,polygon_side):
    angles = np.arctan2(edge_y,edge_x) + np.pi

    ls = np.linspace(0,2*np.pi,
                     num=polygon_side,
                     endpoint=False)
    angle_dists = np.abs(angles[:,np.newaxis] - ls[np.newaxis,:])
    angle_dists = np.where(angle_dists >= np.pi,
                           2 * np.pi - angle_dists,
                           angle_dists)
    min_dist = np.argmin(angle_dists,axis=0)

    bounding_polygon = [edge_y[min_dist],
                        edge_x[min_dist]]

    return bounding_polygon

class SegmentationDataset():
    """Carries segmentation datasets."""

    def __init__(self,
                 hdf5_file,
                 dimensions,
                 rel_keys=['image','mask'],
                 n_anchors=0,
                 transform=None,
                 mode='full',
                 output_scale=[(64,64),(32,32),(16,16)],
                 anchors=[[(1,1),(1.5,1.5)],[(1.5,1.5),(1,1)],[(1,1),(1.5,1.5)]],
                 rotate_record=False,
                 n_polygon_sides=8):
        """

        """
        self.hdf5_file = hdf5_file
        self.hf = h5py.File(self.hdf5_file, 'r')
        self.hf_keys = list(self.hf.keys())
        self.idx_to_keys = {i:x for i,x in enumerate(self.hf_keys)}
        self.dimensions = dimensions
        self.size = (self.dimensions[2]-self.dimensions[0],
                     self.dimensions[3]-self.dimensions[1])

        self.rel_keys = rel_keys
        self.transform = transform
        self.mode = mode
        self.output_scale = [np.array(x) for x in output_scale]
        if isinstance(anchors,list):
            self.anchors = [[np.array(y) for y in x] for x in anchors]
        elif isinstance(anchors,dict):
            self.anchors = {x:[np.array(y) for y in anchors[x]] for x in anchors}
        self.rotate_record = rotate_record
        self.n_polygon_sides = n_polygon_sides

    def __len__(self):
        return len(self.hf)

    def keys(self):
        return list(self.hf.keys())

    def getitem(self,idx):
        def get_rotation_dimensions(h,w,angle):
            theta = np.radians(angle)
            inter_h = np.maximum(
                int(np.ceil(np.sin(theta)*h+np.cos(theta)*w)),
                out_h
                )
            inter_w = np.maximum(
                int(np.ceil(np.cos(theta)*h+np.sin(theta)*w)),
                out_w
            )
            return inter_w,inter_h

        p = 1 # offset for two-step rotation

        if isinstance(idx,int):
            key = self.idx_to_keys[idx]
        else:
            key = idx
        record = self.hf[key]

        in_h,in_w,_ = record['image'].shape
        out_h,out_w = self.size

        if self.rotate_record == True:
            angle = np.random.randint(0,90)
            if angle == 0:
                out_x = np.random.randint(0,in_h - out_h)
                out_y = np.random.randint(0,in_w - out_w)
                record = select_region(
                    record,
                    [out_x,out_y,out_x+out_h,out_y+out_w])

            else:
                inter_h,inter_w = get_rotation_dimensions(out_h,out_w,angle)
                inter_x = np.random.randint(p,in_h - inter_h - p)
                inter_y = np.random.randint(p,in_w - inter_w - p)
                prerotation_dim = (
                    inter_x-p,inter_y-p,
                    inter_x + inter_h + p,inter_y + inter_w + p
                )
                postrotation_dim = (
                    ((inter_h+p) - out_h)//2,
                    ((inter_w+p) - out_w)//2,
                    out_h + ((inter_h+p) - out_h)//2,
                    out_w + ((inter_w+p) - out_w)//2
                )
                record = select_region(record,prerotation_dim)
                record = record_rotation(record,angle)
                record = select_region(record,postrotation_dim)

        else:
            out_x = np.random.randint(0,in_h - out_h)
            out_y = np.random.randint(0,in_w - out_w)
            dimensions = (out_x,out_y,out_x+out_h,out_y+out_w)

            record = select_region(record,dimensions)

        sample = {x:record[x] for x in self.rel_keys}
        sample['image_name'] = key

        if self.transform:
            sample = self.transform(sample)

        return sample

    def intersection_boxes(self,box1,box2):
        return 1 - np.any(
            [(box1[:,0] > box2[:,2]),
             (box2[:,0] > box1[:,2]),
             (box1[:,1] > box2[:,3]),
             (box2[:,1] > box1[:,3])],
            axis=0
        )

    def intersection_area_boxes(self,box1,box2):
        intersects = self.intersection_boxes(box1,box2)
        if np.any(intersects==True):
            x1 = np.where(
                (box1[:,0]<=box2[:,0]) * (box2[:,0]<=box1[:,2]),
                box2[:,0],
                box1[:,0]
            )
            y1 = np.where(
                (box1[:,1]<=box2[:,1]) * (box2[:,1]<=box1[:,3]),
                box2[:,1],
                box1[:,1]
            )
            x2 = np.where(
                box2[:,2]<=box1[:,2],
                box2[:,2],
                box1[:,2]
            )
            y2 = np.where(
                box2[:,3]<=box1[:,3],
                box2[:,3],
                box1[:,3]
            )
            inter_area = np.abs(x2 - x1) * np.abs(y2 - y1)
            union_area = np.add(
                np.abs(box1[:,2]-box1[:,0])*np.abs(box1[:,3]-box1[:,1]),
                np.abs(box2[:,2]-box2[:,0])*np.abs(box2[:,3]-box2[:,1])
            ) - inter_area
            return intersects*(inter_area / union_area)
        else:
            return intersects

    def getitem_object_detection_boxes(self,record):

        sh = [x for x in record['image'].shape[::-1]]

        final_output = {
            'size':{},
            'object':{},
            'class':{},
            'anchors':{}
        }

        if record['centers'].shape[0] > 0:
            bboxes = record['bounding_boxes']
            sizes = np.stack(
                [record['bounding_boxes'][:,2]-record['bounding_boxes'][:,0],
                 record['bounding_boxes'][:,3]-record['bounding_boxes'][:,1]]
            ).T
            centers = np.stack([record['centers'][:,1],
                                record['centers'][:,0]],
                               axis=1)
            c_idxs = np.logical_and(centers[:,0]<(sh[1]-1),
                                    centers[:,1]<(sh[2]-1))

            for scale,anchor_id in zip(self.output_scale,self.anchors):
                anchors = self.anchors[anchor_id]
                outputs_bbox = []
                outputs_object = []
                outputs_class = []
                outputs_anchors = []
                rescaled_sizes = sizes[c_idxs,:] / scale
                rescaled_centers = centers[c_idxs,:] / scale
                rescaled_boxes = np.float32(bboxes[c_idxs,:])
                rescaled_boxes[:,0] = rescaled_boxes[:,0]/scale[0]
                rescaled_boxes[:,1] = rescaled_boxes[:,1]/scale[1]
                rescaled_boxes[:,2] = rescaled_boxes[:,2]/scale[0]
                rescaled_boxes[:,3] = rescaled_boxes[:,3]/scale[1]

                c_x,c_y = [np.floor(rescaled_centers[:,0]).astype(np.int16),
                           np.floor(rescaled_centers[:,1]).astype(np.int16)]

                for anchor in anchors:
                    all_anchor_boxes = np.stack([
                        rescaled_centers[:,0]-anchor[0]/2,
                        rescaled_centers[:,1]-anchor[1]/2,
                        rescaled_centers[:,0]+anchor[0]/2,
                        rescaled_centers[:,1]+anchor[1]/2
                    ],axis=1)
                    sizes_ = np.stack(
                        [rescaled_boxes[:,2]-rescaled_boxes[:,0],
                         rescaled_boxes[:,3]-rescaled_boxes[:,1]]
                    ).T

                    inter = self.intersection_area_boxes(
                        rescaled_boxes,all_anchor_boxes)

                    rescaled_anchored_sizes = np.log(rescaled_sizes / anchor)
                    output_size = sh[1] // scale[0],sh[2] // scale[1]
                    output_object = np.zeros((*output_size,1))
                    output_bbox = np.zeros((*output_size,4))
                    output_class = np.zeros((*output_size,1))
                    output_bbox[c_x,c_y,0] = rescaled_centers[:,0] % 1
                    output_bbox[c_x,c_y,1] = rescaled_centers[:,1] % 1
                    output_bbox[c_x,c_y,2] = rescaled_anchored_sizes[:,0]
                    output_bbox[c_x,c_y,3] = rescaled_anchored_sizes[:,1]
                    output_object[c_x,c_y,0] = inter
                    output_class[c_x,c_y] = 1.
                    output_anchors = np.stack([
                        np.ones(output_size)*anchor[0],
                        np.ones(output_size)*anchor[1]
                    ],axis=-1)

                    outputs_bbox.append(output_bbox)
                    outputs_object.append(output_object)
                    outputs_class.append(output_class)
                    outputs_anchors.append(output_anchors)

                outputs_bbox = np.concatenate(outputs_bbox,axis=-1)
                outputs_object = np.concatenate(outputs_object,axis=-1)
                outputs_class = np.concatenate(outputs_class,axis=-1)
                outputs_anchors = np.concatenate(outputs_anchors,axis=-1)

                scale_id = ','.join([str(x) for x in scale.tolist()])
                final_output['size'][scale_id] = outputs_bbox
                final_output['object'][scale_id] = outputs_object
                final_output['class'][scale_id] = outputs_class
                final_output['anchors'][scale_id] = outputs_anchors

        else:
            for scale,anchor_id in zip(self.output_scale,self.anchors):
                anchors = self.anchors[anchor_id]
                output_size = sh[1] // scale[0],sh[2] // scale[1]
                scale_id = ','.join([str(x) for x in scale.tolist()])
                final_output['size'][scale_id] = np.zeros(
                    (*output_size,len(anchors)*4))
                final_output['object'][scale_id] = np.zeros(
                (*output_size,len(anchors)*1))
                final_output['class'][scale_id] = np.zeros(
                    (*output_size,len(anchors)*1))
                final_output['anchors'][scale_id] = np.concatenate(
                    [np.stack(
                        [
                            np.ones(output_size)*anchor[0],
                            np.ones(output_size)*anchor[1]
                        ],axis=-1
                    ) for anchor in anchors],axis=-1
                )

        final_output['image'] = record['image']

        final_output['image'] = np.transpose(final_output['image'],(2,0,1))
        final_output['size'] = {
            x:np.transpose(final_output['size'][x],(2,0,1))
            for x in final_output['size']}
        final_output['object'] = {
            x:np.transpose(final_output['object'][x],(2,0,1))
            for x in final_output['object']}
        final_output['class'] = {
            x:np.transpose(final_output['class'][x],(2,0,1))
            for x in final_output['class']}
        final_output['anchors'] = {
            x:np.transpose(final_output['anchors'][x],(2,0,1))
            for x in final_output['anchors']}

        return final_output

    def getitem_object_detection_polygons(self,record):

        def get_polygon_distances(centers,polygon_sides,size):
            #size += 1
            max_sizes = np.sqrt((size[:,0])**2 + (size[:,1])**2)/2
            max_angles = np.arctan(size[:,0] / size[:,1])
            cos_ls = np.cos(ls)
            cos_ls_2 = np.cos(np.pi/2 - ls)
            cos_ls[cos_ls < 1e-8] = 1
            cos_ls_2[cos_ls_2 < 1e-8] = 1
            up_dist = 0.5*size[:,0][:,np.newaxis] / cos_ls
            side_dist = 0.5*size[:,1][:,np.newaxis] / cos_ls_2

            sizes = np.where(up_dist <= (max_sizes[:,np.newaxis]),
                             up_dist,
                             side_dist)
            sizes = np.minimum(up_dist,side_dist)
            sizes = np.where(ls <= max_angles[:,np.newaxis],
                             side_dist,
                             up_dist)

            output = []

            for c,ps in zip(centers,polygon_sides):
                D = np.sqrt(np.sum((ps.T-c[np.newaxis,:])**2,axis=1))
                output.append(D)

            output_ = np.array(output)
            output = output_ / sizes
            """
            if np.all(output <= 1.) == False:
                print("FLAG")
                for o,o_,si in zip(output,output_,sizes):
                    for i,(a,c,b) in enumerate(zip(o_,si,o)):
                        if b >= 1.0:
                            print(i,a-c,b)
                        else:
                            print(i)
            """

            output = np.where(output > 1.,1.,output)
            return output

        sh = [x for x in record['image'].shape[::-1]]

        final_output = {
            'size':{},
            'object':{},
            'class':{},
            'polygon':{},
            'anchors':{}
        }

        if record['centers'].shape[0] > 0:
            ls = np.linspace(
                0,2*np.pi,self.n_polygon_sides,
                endpoint=False,dtype=np.float64)[np.newaxis,:] % np.pi
            ls = np.where(
                ls <= (np.pi / 2),ls,np.pi-ls)
            polygon_boxes = np.array([
                record['bounding_polygons'][str(self.n_polygon_sides)][:,0,:].min(axis=1),
                record['bounding_polygons'][str(self.n_polygon_sides)][:,1,:].min(axis=1),
                record['bounding_polygons'][str(self.n_polygon_sides)][:,0,:].max(axis=1),
                record['bounding_polygons'][str(self.n_polygon_sides)][:,1,:].max(axis=1),
            ]).T
            bboxes = polygon_boxes
            sizes = np.stack(
                [np.subtract(bboxes[:,2],bboxes[:,0]),
                 np.subtract(bboxes[:,3],bboxes[:,1])]
            ).T
            centers = record['centers']
            centers = np.stack(
                [np.add(bboxes[:,2],bboxes[:,0])/2,
                 np.add(bboxes[:,3],bboxes[:,1])/2]
            ).T
            c_idxs = np.logical_and(centers[:,0]<(sh[1]-1),
                                    centers[:,1]<(sh[2]-1))

            rescaled_polygon_distances = get_polygon_distances(
                    centers,record['bounding_polygons'][str(self.n_polygon_sides)],sizes
            )
            rescaled_polygon_distances = rescaled_polygon_distances[c_idxs,:]

            for scale,anchors in zip(self.output_scale,self.anchors):
                outputs_bbox = []
                outputs_object = []
                outputs_class = []
                outputs_polygon = []
                outputs_anchors = []
                rescaled_sizes = sizes[c_idxs,:] / scale
                rescaled_centers = centers[c_idxs,:] / scale
                rescaled_boxes = np.float32(bboxes[c_idxs,:])
                rescaled_boxes[:,0] = rescaled_boxes[:,0]/scale[0]
                rescaled_boxes[:,1] = rescaled_boxes[:,1]/scale[1]
                rescaled_boxes[:,2] = rescaled_boxes[:,2]/scale[0]
                rescaled_boxes[:,3] = rescaled_boxes[:,3]/scale[1]

                c_x,c_y = [np.floor(rescaled_centers[:,0]).astype(np.int16),
                           np.floor(rescaled_centers[:,1]).astype(np.int16)]

                for anchor in anchors:
                    all_anchor_boxes = np.stack([
                        rescaled_centers[:,0]-anchor[0]/2,
                        rescaled_centers[:,1]-anchor[1]/2,
                        rescaled_centers[:,0]+anchor[0]/2,
                        rescaled_centers[:,1]+anchor[1]/2
                    ],axis=1)
                    sizes_ = np.stack(
                        [rescaled_boxes[:,2]-rescaled_boxes[:,0],
                         rescaled_boxes[:,3]-rescaled_boxes[:,1]]
                    ).T

                    inter = self.intersection_area_boxes(
                        rescaled_boxes,all_anchor_boxes)

                    rescaled_anchored_sizes = np.log(rescaled_sizes / anchor)
                    #rescaled_anchored_sizes = rescaled_sizes
                    output_size = sh[1] // scale[0],sh[2] // scale[1]
                    output_object = np.zeros((*output_size,1))
                    output_bbox = np.zeros((*output_size,4))
                    output_class = np.zeros((*output_size,1))
                    output_polygon_distances = np.zeros((*output_size,self.n_polygon_sides))
                    output_bbox[c_x,c_y,0] = rescaled_centers[:,0] % 1
                    output_bbox[c_x,c_y,1] = rescaled_centers[:,1] % 1
                    output_bbox[c_x,c_y,2] = rescaled_anchored_sizes[:,0]
                    output_bbox[c_x,c_y,3] = rescaled_anchored_sizes[:,1]
                    output_object[c_x,c_y,0] = inter
                    output_class[c_x,c_y] = 1.
                    output_polygon_distances[c_x,c_y] = rescaled_polygon_distances
                    output_anchors = np.stack([
                        np.ones(output_size)*anchor[0],
                        np.ones(output_size)*anchor[1]
                    ],axis=-1)

                    outputs_bbox.append(output_bbox)
                    outputs_object.append(output_object)
                    outputs_class.append(output_class)
                    outputs_polygon.append(output_polygon_distances)
                    outputs_anchors.append(output_anchors)

                outputs_bbox = np.concatenate(outputs_bbox,axis=-1)
                outputs_object = np.concatenate(outputs_object,axis=-1)
                outputs_class = np.concatenate(outputs_class,axis=-1)
                outputs_polygon = np.concatenate(outputs_polygon,axis=-1)
                outputs_anchors = np.concatenate(outputs_anchors,axis=-1)

                scale_id = ','.join([str(x) for x in scale.tolist()])
                final_output['size'][scale_id] = outputs_bbox
                final_output['object'][scale_id] = outputs_object
                final_output['class'][scale_id] = outputs_class
                final_output['polygon'][scale_id] = outputs_polygon
                final_output['anchors'][scale_id] = outputs_anchors

        else:
            for scale,anchors in zip(self.output_scale,self.anchors):
                output_size = sh[1] // scale[0],sh[2] // scale[1]
                scale_id = ','.join([str(x) for x in scale.tolist()])
                final_output['size'][scale_id] = np.zeros(
                    (*output_size,len(anchors)*4))
                final_output['object'][scale_id] = np.zeros(
                (*output_size,len(anchors)*1))
                final_output['class'][scale_id] = np.zeros(
                    (*output_size,len(anchors)*1))
                final_output['polygon'][scale_id] = np.zeros(
                    (*output_size,self.n_polygon_sides*len(anchors)))
                final_output['anchors'][scale_id] = np.concatenate(
                    [np.stack(
                        [
                            np.ones(output_size)*anchor[0],
                            np.ones(output_size)*anchor[1]
                        ],axis=-1
                    ) for anchor in anchors],axis=-1
                )

        final_output['image'] = record['image']

        final_output['image'] = np.transpose(final_output['image'],(2,0,1))
        final_output['size'] = {
            x:np.transpose(final_output['size'][x],(2,0,1))
            for x in final_output['size']}
        final_output['object'] = {
            x:np.transpose(final_output['object'][x],(2,0,1))
            for x in final_output['object']}
        final_output['polygon'] = {
            x:np.transpose(final_output['polygon'][x],(2,0,1))
            for x in final_output['polygon']}
        final_output['class'] = {
            x:np.transpose(final_output['class'][x],(2,0,1))
            for x in final_output['class']}
        final_output['anchors'] = {
            x:np.transpose(final_output['anchors'][x],(2,0,1))
            for x in final_output['anchors']}

        return final_output

    def getitem_segmentation(self,record):
        return {
            'image':record['image'],
            'mask':record['mask'],
            'weight_map':record['weight_map']
            }

    def __getitem__(self,idx):
        record = self.getitem(idx)
        if self.mode == 'full':
            record = record
        elif self.mode == 'segmentation':
            record = self.getitem_segmentation(record)
        elif self.mode == 'object_detection_boxes':
            record = self.getitem_object_detection_boxes(record)
        elif self.mode == 'object_detection_polygons':
            record = self.getitem_object_detection_polygons(record)
        if isinstance(idx,int):
            key = self.idx_to_keys[idx]
        else:
            key = idx
        record['image_name'] = key
        return record
