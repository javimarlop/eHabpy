import csv
import glob
import os
import sys
import ehab as ehab
import basemodel

class Model(basemodel.BaseModel):
    def execute(self, runparams):
        """Executes the model and writes the output to the CSV file specified in the project params.
        """
        # get parameters qualified by input and output file paths
        params = runparams['parameters']

        # Log the variable values
#         self.logger.fine("random.seed: %s" % params["random.seed"])

        ehab.ehabitat(params["ecoreg_id"], params["shared_dir"], params["output_dir2"])
        
#         # write some variables to output files
#         with open(params["output.filename"], 'w') as f1:
#             f1.write(str(params["test.variable.1"]) + "\n")
# 
