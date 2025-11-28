import os
from PIL import Image

d = "images/minibows2/"

os.makedirs(d, exist_ok=True)

for filename in os.listdir(d):
    if filename.lower().endswith('.pgm'):
        pgm_path = os.path.join(d, filename)
        png_filename = os.path.splitext(filename)[0] + '.png'
        png_path = os.path.join(d, png_filename)
        
        with Image.open(pgm_path) as img:
            img.save(png_path)
            os.remove(pgm_path)