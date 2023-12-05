#!/usr/bin/env python

# import any necessary functions
# argparse is always required
import argparse

# only import these functions when using from template import *
# fill this in with any functions that are not "parse_args()" or "main()"
__all__ = ["my_function"]


# parse command line arguments
def parse_args():
    # initialize an argument parser object
    parser = argparse.ArgumentParser()

    # add arguments one at a time to that argument parser object
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Example required argument accepting a single argument.",
    )
    parser.add_argument(
        "-v",
        "--variable-arg",
        nargs="?",
        default="",
        help="Example optional argument that accepts 0 or 1 arguments.",
    )
    parser.add_argument(
        "-a",
        "--arg-list",
        nargs="+",
        default=[],
        help=(
            "Example optional argument that accepts 1+ space-separated arguments. "
            "By default, returns an empty list."
        ),
    )

    # actually run the argument parsing and accept input from the command line
    args = parser.parse_args()

    # return valid arguments as an object
    return args


# script-specific function that does something useful
def my_function(input_arg: str, variable_arg="", arg_list=None):
    print("This is what you passed to the interpreter:")
    print(f"--input: {input_arg}")
    print(f"--variable-arg: {variable_arg}")
    print(f"--arg-list: {arg_list}")


# where the actual functions you've defined get called
# this part of the script is actually executed
def main():
    # run arg parsing function and return arguments as object
    args = parse_args()

    # access each of the arguments as a variable
    input_file = args.input
    variable_arg = (
        args.variable_arg
    )  # note: dashes in command line arg names are converted to underscores
    arg_list = args.arg_list

    my_function(input_file, variable_arg=variable_arg, arg_list=arg_list)


# this checks whether the script is called from the command line
# if it is, then runs the main() function
if __name__ == "__main__":
    main()
