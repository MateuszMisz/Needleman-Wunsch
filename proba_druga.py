import copy
import sys
from dataclasses import dataclass
from typing import List,Tuple,Optional
import argparse

import numpy as np


@dataclass()
class Sequence:
    name:str
    sequence:str
class MatrixField:
    def __init__(self):
        self.score:Optional[int]=None
        self.previous_field:Optional[List[Tuple[int,int]]]=[]
def print_matrix(matrix):
    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            print(f"{matrix[row][column].score:>{4}}",end='')
        print()

def print_matrix_backtracking(matrix):
    colored_fields=set()
    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            if matrix[row][column].previous_field is not None:
                colored_fields.update(matrix[row][column].previous_field)
    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            if (row,column) in colored_fields:
                print(f"{matrix[row][column].score:>{4}}",end="")
            else:
                print(f"{matrix[row][column].score:>4}",end="")
        print()
def NoneType_aware_max(numbers):
    not_None_numbers=[]
    for number in numbers:
        if number is not None:
            not_None_numbers.append(number)
    return max(not_None_numbers)
def check_parameters():
    if(len(sys.argv)!=2):
        print("Usage: python main.py <input_file>")
def load_from_file(file_path:str):
    with open(file_path, 'r') as file:
        return file.read()
def extract_to_object(fasta_content:str)->List[Sequence]:
    sequences:List=[]
    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            sequences.append(Sequence(name=line,sequence=""))
        else:
            sequences[len(sequences)-1].sequence+=line
    return sequences

def initialize_matrix(list_of_sequences:List[Sequence])->List[List]:
    row_list=[MatrixField() for char in range(len(list_of_sequences[1].sequence)+1)]
    matrix=[copy.deepcopy(row_list) for char in range(len(list_of_sequences[0].sequence)+1)]
    matrix[0][0].score = 0
    return matrix
def count_across(matrix,row,column):
    if row-1<0 or column-1<0:
        return None
    global scores
    global sequences
    if(sequences[0].sequence[row-1]==sequences[1].sequence[column-1]):
        score=scores['match']+matrix[row-1][column-1].score
    else:
        score=scores['mismatch']+matrix[row-1][column-1].score
    return score
def count_from_left(matrix,row,column):
    if column-1<0:
        return None
    global scores
    return scores['gap']+matrix[row][column-1].score
def count_from_top(matrix,row,column):
    if row-1<0:
        return None
    global scores
    return scores['gap']+matrix[row-1][column].score
def fill_cell(matrix,row,column):
    from_across=count_across(matrix,row,column)
    from_left=count_from_left(matrix,row,column)
    from_top=count_from_top(matrix,row,column)
    max_score=NoneType_aware_max([from_across,from_left,from_top])
    matrix[row][column].score=max_score
    matrix[row][column].previous_field=[]
    if max_score==from_across:
        matrix[row][column].previous_field.append((row-1,column-1))
    if max_score==from_left:
        matrix[row][column].previous_field.append((row,column-1))
    if max_score==from_top:
        matrix[row][column].previous_field.append((row-1,column))

def fill_matrix(matrix):
    for row in range(len(matrix)):
        for column in range(len(matrix[0])):
            if row == 0 and column == 0:
                continue
            fill_cell(matrix,row,column)

def backtrack(matrix):
    alignement=[]

    row=len(matrix)-1
    column=len(matrix[0])-1
    while(True):
        tmp_score=[]
        alignement.append((row,column))
        if matrix[row][column].previous_field is None or matrix[row][column].previous_field==[]:
            return alignement
        for field in matrix[row][column].previous_field:
            tmp_score.append(matrix[field[0]][field[1]].score)
        for field in matrix[row][column].previous_field:
            if matrix[field[0]][field[1]].score==NoneType_aware_max(tmp_score):
                row,column=field
def is_going_left(current:Tuple[int,int],next:Tuple[int,int])->bool:
    if current[0]==next[0] and current[1]==next[1]+1:
        return True
    else: return False
def is_going_up(current:Tuple[int,int],next:Tuple[int,int])->bool:
    if current[0]==next[0]+1 and current[1]==next[1]:
        return True
    else:
        return False
def is_going_across(current:Tuple[int,int],next:Tuple[int,int])->bool:
    if current[0]==next[0]+1 and current[1]==next[1]+1:
        return True
    else:
        return False
def backtracking_to_sequences(alignement)->Tuple:
    global sequences
    top_sequence=""
    left_sequence=""
    for i in range(len(alignement)-1):
        if is_going_left(alignement[i],alignement[i+1]):
            left_sequence += "-"
            top_sequence += sequences[1].sequence[alignement[i][1] - 1]

        elif is_going_up(alignement[i],alignement[i+1]):
            top_sequence += "-"
            left_sequence += sequences[0].sequence[alignement[i][0] - 1]
        elif is_going_across(alignement[i],alignement[i+1]):
            left_sequence+=sequences[0].sequence[alignement[i][0]-1]
            top_sequence+=sequences[1].sequence[alignement[i][1]-1]

    return left_sequence[::-1],top_sequence[::-1]

def generate_output(output_sequences:Tuple[str],score:int=0):
    sequences=f"sequence2:\t{output_sequences[0]}\nsequence1:\t{output_sequences[1]}"
    sequences_and_score=f"{sequences}\nscore: {score}"
    return sequences_and_score
def get_output_score(matrix):
    return matrix[len(matrix)-1][len(matrix[0])-1].score
def generate_dev_output(alignement,output_sequences,score:int=0):
    global sequences
    global scores
    return f"""seq1: {sequences[0].sequence}\nseq2: {sequences[1].sequence}
    match: {scores['match']}, mismatch: {scores['mismatch']}, gap: {scores['gap']}
    alignement:{alignement}
    sequence1:\t {output_sequences[0]}
    sequence2:\t {output_sequences[1]}
    score: {score}
    """
def handle_arguments():
    """handles script arguments,
    returns Tuple of input file path and dictionary
    with scores for match,mismatch and gap"""
    argparser = argparse.ArgumentParser(description='needleman_wunsch algorithm')
    argparser.add_argument('input_file', type=str, help='input file(.fasta)')
    argparser.add_argument('--match', type=int, default=1, help='match score')
    argparser.add_argument('--mismatch', type=int, default=-1, help='match score')
    argparser.add_argument('--gap', type=int, default=-2, help='gap score')
    argparser.add_argument('-o','--output', type=str, help='output file',required=True)
    args = argparser.parse_args()
    return args

args=handle_arguments()
file_path=args.input_file
scores={'match':args.match,'mismatch':args.mismatch,'gap':args.gap}
output_path=args.output
fasta_content=load_from_file(file_path)
sequences=extract_to_object(fasta_content)
matrix=initialize_matrix(sequences)
fill_matrix(matrix)
print_matrix_backtracking(matrix)
alignement=backtrack(matrix)
output_sequences=backtracking_to_sequences(alignement)
with open(output_path,'w') as output_file:
    output_file.write(generate_output(output_sequences,get_output_score(matrix)))
print(generate_dev_output(alignement,output_sequences,get_output_score(matrix)))