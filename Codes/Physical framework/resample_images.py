import rasterio as rs
from rasterio.enums import Resampling
import numpy as np
import os


class ResampleImages:
    # Class for resampling images of a given flight
   
    def __init__(self, paths, resample_type):
        # Initializes the paths of the images, the resample type (up or down), resample state and scale factor
        self.paths = paths
        self.resample_type = resample_type
        self.resample_state = self.resample_images_state()
        self.scale_factor = self.get_scale_factor()
   
    def resample_images_state(self):
        # Returns a boolean dictionary with accordance to the resample type:
        # 'downsample' returns True for RGB and DSM; 'upsample' returns True for thermal
        resample_state = np.array([True, True, False])
        if self.resample_type == 'upsample':
            resample_state = ~ resample_state
        return dict(zip(self.paths, resample_state))
       
    def get_scale_factor(self):
        rgb_shape = np.array(rs.open(self.paths[0]).shape)
        ir_shape = np.array(rs.open(self.paths[-1]).shape)
        rgb_ir_ratio = ir_shape / rgb_shape
        if self.resample_type == 'upsample':
            return 1 / rgb_ir_ratio
        return rgb_ir_ratio
   
    def resample_images(self, resamples_dir):
        # Resamples the images of self.paths according to the resample type and returns a list of updated paths
        updated_paths = []
        for path in self.paths:
            if not self.resample_state[path]:
                updated_paths.append(path)
                continue
            image_filename = '/'.join(path.split('/')[4:])[:-4]
            resampled_path = os.path.join(resamples_dir, f'{image_filename}_{self.resample_type}d.tif')
            if os.path.isfile(resampled_path):
                updated_paths.append(resampled_path)
                continue
            with rs.open(path) as image:
                resampled = image.read(out_shape=
                                       (image.count, int(image.height*self.scale_factor[0]),
                                        int(image.width*self.scale_factor[1])),
                                       resampling=Resampling.lanczos, masked=True)
                transform = image.transform*image.transform.scale(image.width/resampled.shape[-1],
                                                                  image.height/resampled.shape[-2])
                profile = image.profile
                profile.update({'height': int(image.height*self.scale_factor[0]),
                                'width': int(image.width*self.scale_factor[1]),
                                'transform': transform})
                updated_paths.append(resampled_path)
                os.makedirs(os.path.dirname(resampled_path), exist_ok=True)
                with rs.open(resampled_path, 'w', **profile) as res_image:
                    res_image.write(resampled)
        return updated_paths
