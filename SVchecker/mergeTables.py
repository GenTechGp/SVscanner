import sys
import pandas as pd
import argparse
import os
import ast

# Takes a list of tab-delimited text files and merges them into one based on a common ID (specified in column number).
# The ID should be unique within each file. Useful for grouping lists of values where each list may not contain every ID.

def mergeOnIndex(cID, files):
    # python3 mergeTab.py testTab/test5.tab testTab/test4.tab -c 0
    dataframes = []

    # Save the first file into the final df
    mainDF = pd.read_csv(files[0], delimiter='\s+', header=None)
    mainDF  = mainDF.astype(str) # convert all values to str
    
    # Read each file into a df and add to list
    for file in files[1:]:
        addDF = pd.read_csv(file, delimiter='\s+', header=None)
        addDF = addDF.astype(str) # convert all values to str
       
        # Check if column index exists in the file - must be less than number of columns
        if cID >= len(addDF.columns):
            sys.stderr.write(f"Error: Could not merge '{file}' - does not contain index {cID}\n")
            continue
        
        dataframes.append(addDF)

    # Merge the other df into the final df using outer join 
    for idx, df in enumerate(dataframes):
        mainDF = pd.merge(mainDF, df, on=cID, how='outer')
    
    return mainDF

def mergeOnDiffIndex(i1, i2, files):
    # python3 mergeTab.py testTab/test8.tab testTab/test9.tab -ui -i1 "[0, 2]" -i2 "[0, 1]"
    # python3 mergeTab.py testTab/test9.tab testTab/test10.tab -ui -i1 3 -i2 0
    # print("Unique index merge")
    mainDF = pd.read_csv(files[0], delimiter='\s+', header=None)
    DF2 = pd.read_csv(files[1], delimiter='\s+', header=None)
    
    # Argument is a list
    if i1.startswith('[') and i1.endswith(']'):
        i1 = ast.literal_eval(i1)
    else: 
        i1 = int(i1)
    if i2.startswith('[') and i2.endswith(']'):
        i2 = ast.literal_eval(i2)
    else:
        i2 = int(i2)


    mainDF = pd.merge(mainDF, DF2, how='left', left_on=i1, right_on=i2)

    return mainDF

    #     # Merge the depth values for the start positions
    # merged_start_df = pd.merge(reads_df, depth_df, how='left', left_on=['chromosome', 'start'], right_on=['chromosome', 'position'])
    # merged_start_df = merged_start_df.drop(columns=['position'])

    # # Merge the depth values for the end positions
    # merged_end_df = pd.merge(merged_start_df, depth_df, how='left', left_on=['chromosome', 'end'], right_on=['chromosome', 'position'])
    # merged_df = merged_end_df.drop(columns=['position'])
    

def mergeOnName(nID, files):
    # python3 mergeTab.py testTab/test1.tab testTab/test2.tab -n ID
    dataframes = []

    # Save the first file into the final df
    mainDF = pd.read_csv(files[0], delimiter='\s+')
    mainDF  = mainDF.astype(str) # convert all values to str

    # Read each file into a df and add to list
    for file in files[1:]:
        addDF = pd.read_csv(file, delimiter='\s+')
        addDF = addDF.astype(str) # convert all values to str

        # Check if column header exists in the file 
        if nID not in addDF.columns:
            sys.stderr.write(f"Error: Could not merge '{file}' - does not contain header name {nID}\n")
            continue

        dataframes.append(addDF)

    # Merge the other df into the final df using outer join 
    for idx, df in enumerate(dataframes):
        mainDF = pd.merge(mainDF, df, on=nID, how='outer')

    return mainDF

def mergeOnUnique(n1, n2, files):
    # python3 mergeTab.py testTab/test6.tab testTab/test7.tab -u -n1 "['chr', 'pos']" -n2 "['chr', 'start']"
    # python3 mergeTab.py testTab/test6.tab testTab/test7.tab -u -n1 pos -n2 start <-- Fix drop duplicates
    # Check Two Files
    mainDF = pd.read_csv(files[0], delimiter='\s+')
    DF2 = pd.read_csv(files[1], delimiter='\s+')
    
    dropColumns = []

    # Argument is a list
    if n1.startswith('[') and n1.endswith(']'):
        n1 = ast.literal_eval(n1)
    if n2.startswith('[') and n2.endswith(']'):
        n2 = ast.literal_eval(n2)
        dropColumns = n2[1:]


    mainDF = pd.merge(mainDF, DF2, how='outer', left_on=n1, right_on=n2)
    mainDF = mainDF.drop(columns=dropColumns)

    return mainDF




if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="Merge Tables",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help='tabfile1.tab tabfile2.tab')

    parser.add_argument("-o",
        help="output file name [stdout]")
    parser.add_argument("-c", "--columnID",
        help="column number of common ID (int)")
    parser.add_argument("-head", "--head", action="store_true", 
        help="First row contains a header")
    parser.add_argument("-n", "--headerID",
        help="column header of common ID (string)")
    parser.add_argument("-b", "--blank",
        help="Replace blank fields with specified character", default=0)

    # Create an argument group for the -u flag and its options
    uniqueGroup = parser.add_argument_group('unique Header')

    # Add the -u flag
    uniqueGroup.add_argument('-u', '--uniqueHeader', action='store_true', help='Merge columns with different header names')
    # Add options for column names of the first and second files
    uniqueGroup.add_argument('-n1', '--headerID1', help='Column name for the first file', required='-u' in sys.argv)
    uniqueGroup.add_argument('-n2', '--headerID2', help='Column name for the second file', required='-u' in sys.argv)

    uniqueGroup.add_argument('-ui', '--uniqueIndex', action='store_true', help='Merge columns with different header indicies')
    # Add options for column names of the first and second files
    uniqueGroup.add_argument('-i1', '--columnID1', help='Column name for the first file', required='-ui' in sys.argv)
    uniqueGroup.add_argument('-i2', '--columnID2', help='Column name for the second file', required='-ui' in sys.argv)

    args = parser.parse_args()
    
    isHeader = False
    # mainDF = None

    ##### ERROR CHECKING #####
    for file in args.files:
        # Check at least 2 files 
        if len(args.files) < 2:
            sys.stderr.write('Error: At least two files are required.\n')
            sys.exit(1)
        #  Check filenames are valid
        if not os.path.isfile(file):
            sys.stderr.write(f'Error: File "{file}" does not exist.\n')
            sys.exit(1)

    ##### MERGE TABLES #####
    # Merge on common ID column number
    if args.columnID:
        # print(f'Merging on{args.columnID}')
        mainDF = mergeOnIndex(int(args.columnID), args.files)
    # Merge on common ID column name
    elif args.headerID:    
        mainDF = mergeOnName(args.headerID, args.files)
        isHeader = True
    
    # Merge on two different ID (only two arguments passed)
    # https://saturncloud.io/blog/how-to-perform-a-pandas-join-on-columns-with-different-names/#:~:text=To%20join%20two%20DataFrames%20on,()%20or%20join()%20function.
    elif args.uniqueHeader:
        mainDF = mergeOnUnique(args.headerID1, args.headerID2, args.files)
    elif args.uniqueIndex:
        mainDF = mergeOnDiffIndex(args.columnID1, args.columnID2, args.files)
       
    # # Replace the NaN values with 0 
    mainDF = mainDF.convert_dtypes()

    if args.o:
        mainDF.to_csv(args.o, sep='\t', index=False, header = isHeader)
    else:
        print(mainDF.to_string(index=False, header = isHeader))

    

