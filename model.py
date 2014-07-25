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
        self.logger.fine("Output directory: %s" % params["output_dir"])
        od = params["output_dir"][:-1]
        self.logger.fine("New output directory: %s" % od)
        
        d = os.path.join(os.path.sep, od,'results')
        
        try:
            os.makedirs(d)
        except OSError:
            if os.path.isdir(d):
                self.logger.fine("Can't create directory: %s - a directory with this name already exists. It will be used for the results of the analysis." % d)  
            else: 
                self.logger.fine("Can't create directory: %s - halting the analysis as there is nowhere to write results." % d)  
        
        ehab.ehabitat(params["ecoreg_id"], params["shared_dir"], od)
        
#         # write some variables to output files
#         with open(params["output.filename"], 'w') as f1:
#             f1.write(str(params["test.variable.1"]) + "\n")
# 
