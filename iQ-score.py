import os
import multiprocessing
from tqdm import tqdm
import get_good_inter_pae
import pandas as pd
import csv
import traceback
import argparse



#Path_ccp4 = "/opt/xtal/ccp4-9"
parser = argparse.ArgumentParser()
parser.add_argument("--N_CPU", help="Number of CPUs available for computation ", required=True, type=int)
parser.add_argument("--Path_ccp4", help="Path of ccp4", required=True, type=str)
parser.add_argument("--model_dir", help="Path of interactions directory", required=True, type=str)
parser.add_argument("--msa_dir", help="Path of all MSA", required=True, type=str)
args = parser.parse_args()

def run_scoring(args) :
    """
    Run get_good_inter_pae script.

    Parameters:
    ----------
    args : set

    Returns:
    ----------
    result : string
    """
    interaction, Path_ccp4, seq_no_SP = args
    try:
        result = get_good_inter_pae.main(interaction, 10, 2, Path_ccp4, seq_no_SP) #normal PAE is 10
        return result
    except Exception as e:
        pid = os.getpid()
        print(f"ERROR in worker PID={pid}")
        print(f"Interaction: {interaction}")
        traceback.print_exc()
        raise




def Score_interaction (msa_dir,model_dir,Path_ccp4,N_CPU) :
    """
    Generate scores for all interactions and set a list of possible prey based on the scores.

    Parameters:
    ----------
    file : object of class File_proteins
    model_dir : string
    Informations_dict : dictionary
    Interaction : string

    """
    seq_no_SP = dict()
    ppi_list = list()

    for dir in os.listdir(model_dir) :
        if "_and_" in dir :
            ppi_list.append(model_dir + "/" + dir)
            prot1 = dir.split("_and_")[0]
            with open(f"{msa_dir}/{prot1}.a3m", "r") as a3m1 :
                i = 0
                for line in a3m1 :
                    if i == 1 :
                        seq_no_SP[prot1] = len(line.strip())
                        break
                    if line[0] == ">" :
                        i = 1
            prot2 = dir.split("_and_")[1]
            with open(f"{msa_dir}/{prot2}.a3m", "r") as a3m2 :
                i = 0
                for line in a3m2 :
                    if i == 1 :
                        seq_no_SP[prot2] = len(line.strip())
                        break
                    if line[0] == ">" :
                        i = 1
            Interaction = "PPI_int"
    results = []
    with multiprocessing.Pool(N_CPU) as pool :
         tasks = [(ppi, Path_ccp4, seq_no_SP) for ppi in ppi_list]
         results_iter = pool.imap_unordered(run_scoring, tasks)
         for df in tqdm(results_iter, total=len(ppi_list), desc="Scoring interactions") :
             if df is not None and not df.empty :
                 results.append(df)
         pool.close()
         pool.join()
    if results :
         merged_df = pd.concat(results, ignore_index=True)
         merged_df.to_csv(os.path.join(f"./", "predictions_with_good_interpae.csv"), index=False)

    with open(f"./predictions_with_good_interpae.csv", "r") as result_file :
        reader = csv.DictReader(result_file)

        if Interaction == "PPI_int" : 
            all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
            for row in reader :
                job = row['jobs'].split("ranked")[0]
                if '_and_' in job  : 
                    if row['pi_score'] == 'No interface detected' :
                        iQ_score = float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30 #pi_score don't detect interface so it's set on -2.63
                        line =f'{row["jobs"]},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                    else :
                        iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                        line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'

                    all_lines = all_lines + line
            for ppi in ppi_list :
                line = f"{ppi.split('/')[-1]}_ranked_0,0,0,0,0\n"
                all_lines = all_lines + line
   
    if len(all_lines.strip("\n")) > 1 : #if all_lines is not empty
            with open(f"./iQ-score_results.csv", "w") as file2 :
                file2.write(all_lines)

Score_interaction(args.msa_dir,args.model_dir,args.Path_ccp4, args.N_CPU)
