import urllib
import json

with open('constants.json') as data_file:    
    const = json.load(data_file)

data = urllib.URLopener()
data.retrieve(const["fileURL"], const["outPath"]+const["fileName"])
