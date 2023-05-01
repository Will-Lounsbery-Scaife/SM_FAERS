import os
import time
   
startTime = time.time()

os.system("python scripts/post_analysis.py")

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))