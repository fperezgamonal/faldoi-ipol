
from os import listdir
from os.path import isfile, join


def list_images_dataset(dataset, random_subset):
	if random_subset:
		file_images_random = '../scripts_python/Random_images.txt'
		with open(file_images_random, 'r') as file:
		# read a list of lines into data
			images_dir = file.readlines()
		for i in range(len(images_dir)):
			images_dir[i] = images_dir[i][:-1]
	else:
		images_dir = []
		if dataset == 'sintel':
			directories = ['../data/Sintel_final', '../data/Sintel_clean']
		elif dataset == 'middlebury':
			directories = ['../data/Middlebury']
	
		for directory_images in directories:
		
			folders_sequences = listdir(directory_images)
		
			for i in range(len(folders_sequences)):

				directory_sequence = join(directory_images, folders_sequences[i])
				images = sorted([f for f in listdir(directory_sequence) if isfile(join(directory_sequence, f)) and f.endswith(".txt")])

				images_dir = images_dir + [join(directory_sequence, f) for f in images]#[1:-1]]
			
	return images_dir
