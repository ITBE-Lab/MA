from MA import *
from MSV import *

def compute_accuracy_recall(dataset_name, blurs, gt_ids, run_ids, out_prefix):
    db_conn = DbConn(dataset_name)
    count_calls_from_db = SvCallsFromDb(db_conn)
    for blur in blurs:
        for gt_id in gt_ids:
            for run_id in run_ids[gt_id]:
                print("computing", blur, gt_id, run_id)
                stats, gt_total = count_calls_from_db.count(run_id, gt_id, blur)
                with open(out_prefix + dataset_name + "-" + str(run_id) + "-" + str(blur) + ".tsv", "w") as out_file:
                    out_file.write("//|ground truth| = " + str(gt_total) + "\n")
                    out_file.write("//#supporting reads\t#true positives\t#num entries\trecall\taccuracy\n")
                    for x, num_calls, num_tp in stats:
                        if num_calls > 0:
                            out_file.write(str(x) + "\t" + str(num_tp) + "\t" + str(num_calls) +
                                            "\t" + str(num_tp/gt_total) + "\t" + str(num_tp/num_calls) + "\n")