import urllib
import json

with open('constants.json') as data_file:    
    const = json.load(data_file)

for n in range(len(const["fileURL"])):
    data = urllib.URLopener()
    data.retrieve(const["fileURL"][n], const["outPath"]+const["fileName"][n])

