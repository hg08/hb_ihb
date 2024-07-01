#!/usr/bin/env python
import sys 

if len(sys.argv) != 2:
    print("Ooops. Usage:\n" + 'python ' + sys.argv[0] + ' 128w-pos-1.xyz')
    exit()

# Input file
f1 = open(sys.argv[1],'r')

# Smaple every n frames
ns = 80

# Output file
out_name = sys.argv[1].split('.xyz')[0] + '_sampled.xyz'
f2 = open(out_name, 'w')


# File as a variable to get the atom numbers
f3 = open(sys.argv[1], 'r')
num_atom = int(f3.readline())
f3.close()
print('Number of atoms:{}'.format(num_atom))


def read_one_frame(f1, num_atom):
    '''
    Read one frame from the input file
        Input:
            f1: file object
            num_atom: number of atoms 
        Output:
            i: frame number
            frame: a list of lines in one frame
    '''
    try:
        frame = []
        for _ in range(num_atom + 2):
            line = f1.readline()
            frame.append(line)
        i = frame[1].split(',')[0].split('=')[-1] 
        return int(i), frame
    except:
        print('End of file')
        return None, None

# Read the first frame
i, frame = read_one_frame(f1, num_atom)
while i is not None and i < 80000:
    if i % ns == 0:
        for k, line in enumerate(frame):
            if k == 1:
                field1 = line.split(',')[0]
                line = ' i =         ' + line[len(field1):]
                new_i = str(i // ns)
                line = line[0:len(field1)-len(new_i)] + new_i + line[len(field1):]
            f2.write(line)
    i, frame = read_one_frame(f1, num_atom)

f1.close()
f2.close()

print('{} saved.'.format(out_name))
