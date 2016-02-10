#! /usr/bin/env python
__author__ = 'heroico'

import logging
import os
import re
import KeyedDataSet
import WeightDBUtilities
import GWASUtilities
import MethodGuessing
import Utilities
import Logging

class GetBetas(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.gwas_folder = args.gwas_folder
        self.output_folder = args.output_folder
        self.compressed = args.compressed
        self.args = args
        self.gwas_regexp = None
        if args.gwas_file_pattern:
            self.gwas_regexp = re.compile(args.gwas_file_pattern)

    def run(self):
        logging.info("Loading weight model")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db_path)

        names = Utilities.contentsWithRegexpFromFolder(self.gwas_folder, self.gwas_regexp)

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        for name in names:
            self.buildBetas(weight_db_logic,name)

    def buildBetas(self, weight_db_logic, name):
        output_path = os.path.join(self.output_folder, name)
        if os.path.exists(output_path):
            logging.info("%s already exists, delete it if you want it to be done again", output_path)
            return

        logging.info("Building beta for %s and %s", name, self.weight_db_path)
        input_path = os.path.join(self.gwas_folder, name)
        file_format = GWASUtilities.GWASFileFormat.fileFormatFromArgs(input_path, self.args)

        scheme = MethodGuessing.chooseGWASProcessingScheme(file_format, weight_db_logic, self.args, input_path)

        callback = GWASUtilities.GWASWeightDBFilteredBetaLineCollector(file_format, scheme, weight_db_logic)
        dosage_loader = GWASUtilities.GWASDosageFileLoader(input_path, self.compressed, self.args.separator, callback)
        result_sets = dosage_loader.load()
        KeyedDataSet.KeyedDataSetFileUtilities.saveSetsToCompressedFile(output_path, result_sets, "rsid")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build betas from GWAS data.')

    parser.add_argument("--weight_db_path",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    parser.add_argument("--output_folder",
                        help="name of folder to put results in",
                        default="intermediate/beta")

    parser.add_argument("--scheme",
                        help="Type of beta data preprocessing, optional. Options are: "
                        "'beta' (provide (beta or OR));"
                        "'beta_se' (provide (beta or OR) and standard error); "
                        "'beta_se_to_z' (provide (beta or OR) and standard error), and Zscore of beta will be output;"
                        "'z' (provide zscore of beta),"
                        " 'beta_sign_p' (sign of beta, and pvalue); beta_p (beta and pvalue)",
                        default=None)

    parser.add_argument("--or_column",
                    help="Name of column containing Odd Ratios in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--pvalue_column",
                    help="Name of column containing p-value in input files.",
                    default=None)

    parser.add_argument("--beta_sign_column",
                    help="Name of column containing sign of beta in input files.",
                    default=None)

    parser.add_argument("--beta_column",
                    help="Name of column containing betas in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--se_column",
                    help="Name of column containing standard error in input file.",
                    default=None)

    parser.add_argument("--beta_zscore_column",
                    help="Name of column containing beta's zscore in input file.",
                    default=None)

    parser.add_argument("--frequency_column",
                    help="Name of column containing frequency in input file",
                    default=None)

    parser.add_argument("--a1_column",
                    help="Name of column containing allele 1 in input file (reference allele, following PrediXcan format, and plink --dosage format philosophy)",
                    default="A1")

    parser.add_argument("--a2_column",
                    help="Name of column containing allele 2 in input file (dosage/effect allele)",
                    default="A2")

    parser.add_argument("--snp_column",
                    help="Name of column containing snp in input file",
                    default="SNP")

    parser.add_argument("--compressed",
                    help="Wether input files are gzip compressed file",
                    action="store_true",
                    default=False)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = GetBetas(args)
    if args.throw:
        work.run()
    else:
        try:
            work.run()
        except NameError as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)
        except Exception as e:
            print e
            exit(1)
