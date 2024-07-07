#!/usr/bin/env python
import sys 

if len(sys.argv) != 4:
    print("Ooops. Usage:\n" + 'python ' + sys.argv[0] + ' <original_traj>' + ' <start_frame>' + ' <sampled_traj>')
    exit()

# Original trajectory file
f1 = open(sys.argv[1],'r')
# Start frame
start_frame = int(sys.argv[2]) 
# Output sampled file
out_name = sys.argv[3]
f2 = open(out_name, 'w')


# Smaple every n frames
ns = 80

# Sub-trajectory length in ps
sub_length_in_ps = 15.0 
sub_length_in_fs = sub_length_in_ps * 1000  
dt = 0.5  # fs 
numFrames = int(sub_length_in_fs / dt)  
end_frame = start_frame + numFrames

# File as a variable to get the atom numbers
f3 = open(sys.argv[1], 'r')
num_atom = int(f3.readline())
f3.close()
print('Number of atoms:{}'.format(num_atom))

totalFrames = 80000

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

for j in range(totalFrames):
    i, frame = read_one_frame(f1, num_atom)
    if i is not None and i >= start_frame and i < end_frame:
        if i % ns == 0:
            for k, line in enumerate(frame):
                if k == 1:
                    field1 = line.split(',')[0]
                    line = ' i =         ' + line[len(field1):]
                    new_i = str((i-start_frame) // ns)
                    line = line[0:len(field1)-len(new_i)] + new_i + line[len(field1):]
                f2.write(line)

f1.close()
f2.close()

print('{} saved.'.format(out_name))
