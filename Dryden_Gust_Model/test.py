import csv


# wide range noise generation
t_w = []
noise1 = []
noise2 = []
noise3 = []
with open('white_noise.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',')
    for row in csv_reader:
        t_w.append(float(row[0]))
        noise1.append(float(row[1]))
        noise2.append(float(row[2]))
        noise3.append(float(row[3]))


def main():
    print(noise1)
main()