import copy
import sys
from dataclasses import dataclass
from typing import List,Tuple,Optional
from colorama import Fore, Back, Style, init


import numpy as np

scores={
    'match':1,
    'mismatch':-1,
    'gap':-2
}
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
                print(f"{Fore.RED}{matrix[row][column].score:>{4}}",end="")
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

def backtrack(matrix:List[List[MatrixField]],row:int=-1,column:int=-1,list_of_aligments=[[]],sequence_index=0)->List[List[Tuple[int,int]]]:
    if(row==-1 or column==-1):
        row=len(matrix)-1
        column=len(matrix[0])-1
    list_of_aligments[sequence_index].append((row,column))
    cnt=0
    index=len(list_of_aligments[sequence_index])
    for field_coordinates in matrix[row][column].previous_field:
        if(cnt>0):
            list_of_aligments.append(copy.deepcopy(list_of_aligments[sequence_index][0:sequence_index+1]))
        backtrack(matrix,field_coordinates[0],field_coordinates[1],list_of_aligments,len(list_of_aligments)-1)
        cnt+=1
    return list_of_aligments




#check_parameters()
#file_path=sys.argv[1]

file_path='test.fasta'
fasta_content=load_from_file(file_path)
sequences=extract_to_object(fasta_content)
matrix=initialize_matrix(sequences)
fill_matrix(matrix)
print_matrix_backtracking(matrix)
list_of_aligments=backtrack(matrix)
print(list_of_aligments)
