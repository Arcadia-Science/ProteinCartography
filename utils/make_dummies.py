import argparse
import subprocess
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, nargs = '+')
    parser.add_argument("-o", "--output", required = True)
    parser.add_argument("-M", "--maximum", default = '0')
    args = parser.parse_args()
    return args

def make_dummies(input_file: str, output_dir: str, maximum: int):
    if output_dir[-1] != '/':
        output_dir = output_dir + '/'
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    input_text = []
    with open(input_file[0], 'r') as f:
        input_text = f.read().splitlines()
    
    if maximum == 0 or maximum > len(input_text):
        pass
    else:
        input_text = input_text[:maximum]
    
    for acc in input_text:
        output_path = output_dir + acc + '.txt'
        subprocess.run(['touch', output_path])

def main():
    args = parse_args()
    input_file = args.input
    output_dir = args.output
    maximum = int(args.maximum)
    make_dummies(input_file, output_dir, maximum)
    
if __name__ == '__main__':
    main()
    