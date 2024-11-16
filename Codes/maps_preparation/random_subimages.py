import rasterio as rs
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import time
import random


class RandomSubimages:
    # Class for creating random square subimages from an image
    
    def __init__(self, path):
        # Initializes the image and the set of coordinates that correspond to the subimages
        self.path = path
        self.image = self.path_to_image()
        self.coordinates = []
    
    def path_to_image(self):
        # Gets the image out of self.path
        with rs.open(self.path) as image:
            return image.read()
    
    def is_expandable(self, coordinate, size):
        # Checks if a given coordinate is expandable to a square of size self.size that is contained in the image
        return self.image[3, coordinate[0]:coordinate[0] + size, coordinate[1]:coordinate[1] + size].all()
    
    def is_congruential(self, coordinate, size, cong_percent):
        # Checks if a given coordinate is congruent with the coordinates in self.coordinates with percentage of less than cong_percent
        for c in self.coordinates:
            if (size - abs(coordinate[0] - c[0]))*(size - abs(coordinate[1] - c[1]))/size**2 > cong_percent:
                return False
        return True
    
    def resize(self, size, image):
        # Returns minimum of int(image.shape / 100) * 100 in case size is bigger; Otherwise returns size
        min_image_size = min(image.shape)
        if min_image_size < size:
            print(f'Subimages resized to {int(min_image_size / 100) * 100}')
            return int(min_image_size / 100) * 100
        return size
    
    def create_subimages(self, num, size, cong_percent=0.2, plot=False):
        # Creates a given number of subimages with a given size with percentage of congruence of less than cong_percent
        with rs.open(self.path) as image:
            time_0 = time.time()
            current_cong_percent = cong_percent
            size = self.resize(size, image)
            while len(self.coordinates) < num:
                coordinate = list(map(lambda x: random.randint(0, x), [image.height - size, image.width - size]))
                if self.is_expandable(coordinate, size) and self.is_congruential(coordinate, size, current_cong_percent):
                    self.coordinates.append(coordinate)
                if time.time() - time_0 > 10:
                    print(current_cong_percent)
                    current_cong_percent += 0.05
                    time_0 = time.time()
        if plot:
            self.plot_subimages(size)
        return self.coordinates
    
    def plot_subimages(self, size):
        # Plots the subimages from given coordinates and size on the image
        fig, ax = plt.subplots(figsize=(8,8))
        ax.imshow(self.image.transpose((1, -1, 0)))
        squares = [Rectangle(coordinate[::-1], size, size, edgecolor='r', facecolor='none') for coordinate in self.coordinates]
        for square in squares:
            ax.add_patch(square)
        plt.show()
