import ROOT
import root_pandas
import pandas
import glob

DEFAULT_BRANCH_LIST = ['__experiment__', '__run__', '__event__']


def get_branch_list(paths, branch_list):
    df = root_pandas.read_root(paths, columns=branch_list)
    df_string = df.to_string(index=False, header=False, sparsify=False, max_rows=10, col_space=1)
    df_string = df_string.replace("  ",",").replace(" ","") 
    print(df_string)


def writeEventList(df, outfile, columns=DEFAULT_BRANCH_LIST):
    df.to_csv(outfile, index=False, columns=columns)
    print(f"# output written in {outfile}")

def isInDF(df, **kwargs):
    for k in kwargs:
        if not k in df.columns:
            raise KeyError(f"kwarg key: {k} was not found in df.columns={df.columns}")
    query_string = " and ".join(f"{k}=={v}" for k,v in kwargs.items())
    ret = df.query(query_string)
    return len(ret)

def getValsFromChain(chain, branch_list=DEFAULT_BRANCH_LIST):
    return {k:getattr(chain,k) for k in branch_list}




def eListSelector(chain, elist_df, *args, **kwargs):
    vals = getValsFromChain(chain)
    return isInDF(elist_df, **vals)



if __name__ == "__main__":


    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=[], nargs="+")
    parser.add_argument("--output", default="test_event_run_info.txt")
    #parser.add_argument("--tree_name", default="tau3x1")
    #parser.add_argument("--cut", default='')
    #parser.add_argument("--n", type=int, default=-1)
    parser.add_argument("--columns", default=DEFAULT_BRANCH_LIST, nargs="+", )
    args = parser.parse_args()

    
    paths   = args.input
    if len(paths) == 1 and "*" in paths[0]:
        paths = glob.glob(paths[0])

    columns = args.columns
    
    df = root_pandas.read_root(paths, columns=columns)
    writeEventList(df, args.output, columns=columns)


