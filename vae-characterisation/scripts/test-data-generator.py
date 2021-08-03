from data_generator import DataGenerator

dg = DataGenerator("datasets/wbc.h5")

for d in dg.iterate():
    print(d.shape)
