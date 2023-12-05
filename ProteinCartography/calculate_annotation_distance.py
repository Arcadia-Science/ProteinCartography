#!/usr/bin/env python
import argparse

import nltk
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel

# downloads a pretrained punkt word tokenizer for use
# by nltk.tokenize.word_tokenize
# see https://www.nltk.org/api/nltk.tokenize.punkt.html
nltk.download("punkt")

# only import these functions when using import *
__all__ = ["calculate_annotation_distance"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Path of input features file.")
    parser.add_argument("-o", "--output", required=True, help="Path of output distances file.")
    parser.add_argument(
        "-d", "--id-col", required=True, help="ID column to use for unique entries."
    )
    parser.add_argument("-a", "--annot-col", required=True, help="Column containing annotations.")
    parser.add_argument(
        "-f",
        "--filter-col",
        default=[],
        nargs="?",
        help="Columns to use for filtering.",
    )
    args = parser.parse_args()

    return args


def calculate_annotation_distance(
    input_file: str,
    id_col: str,
    annot_col: str,
    output_file=None,
    sep="\t",
    filter_cols=None,
    dropna=True,
    tokenizer="nltk_word_tokenize",
    ignored_tokens=("(", ")", ";", "."),
    round_digits=8,
    ngram_range=(1, 3),
    min_df=0.0,
):
    """
    Takes a TSV file containing annotations and calculates the all-v-all cosine similarity
    between them.
    The output is a TSV file containing the cosine similarity between each pair of annotations.

    Args:
        input_file (str): path to the input file
        id_col (str): name of the column containing the IDs
        annot_col (str): name of the column containing the annotations
        output_file (str): path to the output file
        sep (str): separator used in the input file
        filter_cols (list): list of columns to be used as filters based on missing values
        dropna (bool): whether to drop rows with missing values in the filter columns
        tokenizer (str): tokenizer to use for tokenizing the annotations
        ignored_tokens (tuple): list of tokens to be ignored
        round_digits (int): number of digits to round the cosine similarity to
        ngram_range (tuple): range of n-grams to be used for the TF-IDF vectorizer
        min_df (float): minimum document frequency for the TF-IDF vectorizer
    """

    if filter_cols is None:
        filter_cols = []

    # Read the input file
    input_df = pd.read_csv(input_file, sep=sep)

    # Select the columns to be used as keys
    annotations = input_df[[id_col, annot_col] + filter_cols]

    # Drop rows with missing values in the filter columns
    if dropna:
        annotations.dropna(subset=filter_cols, inplace=True)

    # Tokenize the annotations
    if tokenizer == "nltk_word_tokenize":
        tokenize = nltk.tokenize.word_tokenize
    elif tokenizer == "nltk_wordpunct_tokenize":
        tokenize = nltk.tokenize.wordpunct_tokenize
    else:
        raise ValueError(
            f"Tokenizer '{tokenizer}' not recognized. "
            "Please use 'nltk_word_tokenize' or 'nltk_wordpunct_tokenize'."
        )
    annotations[annot_col] = annotations[annot_col].apply(tokenize)

    # Remove sanitized words
    annotations[annot_col] = annotations[annot_col].apply(
        lambda x: [word for word in x if word not in ignored_tokens]
    )

    # Create final sanitized string for each annotation
    annotations[annot_col] = annotations[annot_col].apply(lambda x: " ".join(x))

    # Initialize the TF-IDF vectorizer
    # min_df=0.0 ensures that all words are included
    # ngram_range=(1, 3) includes unigrams, bigrams, and trigrams
    tf = TfidfVectorizer(analyzer="word", ngram_range=ngram_range, min_df=min_df)

    # Fit the vectorizer to the sanitized annotations
    tfidf = tf.fit_transform(annotations[annot_col])

    # Calculate the cosine similarity between the sanitized annotations
    # (note: the linear kernel is equivalent to the cosine similarity
    # because the TF-IDF vectors are normalized, and it is faster to calculate)
    tfidf_cosine = linear_kernel(tfidf)
    tfidf_df = pd.DataFrame(tfidf_cosine, columns=annotations[id_col], index=annotations[id_col])
    tfidf_df = tfidf_df.round(round_digits)

    if output_file is not None:
        tfidf_df.to_csv(output_file, sep=sep)

    return tfidf_df


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    id_col = args.id_col
    annot_col = args.annot_col
    filter_cols = args.filter_col

    if isinstance(filter_cols, str):
        filter_cols = [filter_cols]

    calculate_annotation_distance(
        input_file=input_file,
        id_col=id_col,
        annot_col=annot_col,
        output_file=output_file,
        filter_cols=filter_cols,
    )


if __name__ == "__main__":
    main()
