# This is a Python script for retrieving data from the NASA/IPAC Infrared Science Archive (IRSA).
# In this particular code, I'm uploading a list of sources, and opening a new browser tab to the search results from the WISE All-Sky Survey.
# Data on how to use the Gator API can be found here: https://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html

# Import libs
import requests #used to make API request
import xml.etree.ElementTree as ET #used to parse XML to return result data URL
import webbrowser #used to open a new browser tab to result data

# I found that the Gator API only appears to return short form data, and if you want 
# any additional columns, you have to specify those, and to specify those, you have to 
# run a query of the Catalog Data Dictionary Form Service (CatDD) to see a list of column options.
# I ran query once and parsed it into a list object called longcols.
# longcols returns all columns when passed as the selcols query parameter:
longcols=["designation","ra","dec","sigra","sigdec","sigradec","glon","glat",\
	"elon","elat","wx","wy","cntr","source_id","coadd_id","src","w1mpro",\
	"w1sigmpro","w1snr","w1rchi2","w2mpro","w2sigmpro","w2snr","w2rchi2",\
	"w3mpro","w3sigmpro","w3snr","w3rchi2","w4mpro","w4sigmpro","w4snr",\
	"w4rchi2","rchi2","nb","na","w1sat","w2sat","w3sat","w4sat","satnum",\
	"ra_pm","dec_pm","sigra_pm","sigdec_pm","sigradec_pm","pmra","sigpmra",\
	"pmdec","sigpmdec","w1rchi2_pm","w2rchi2_pm","w3rchi2_pm","w4rchi2_pm",\
	"rchi2_pm","pmcode","cc_flags","rel","ext_flg","var_flg","ph_qual","det_bit",\
	"moon_lev","w1nm","w1m","w2nm","w2m","w3nm","w3m","w4nm","w4m","w1cov","w2cov",\
	"w3cov","w4cov","w1cc_map","w1cc_map_str","w2cc_map","w2cc_map_str","w3cc_map",\
	"w3cc_map_str","w4cc_map","w4cc_map_str","best_use_cntr","ngrp","w1flux","w1sigflux",\
	"w1sky","w1sigsk","w1conf","w2flux","w2sigflux","w2sky","w2sigsk","w2conf",\
	"w3flux","w3sigflux","w3sky","w3sigsk","w3conf","w4flux","w4sigflux","w4sky",\
	"w4sigsk","w4conf","w1mag","w1sigm","w1flg","w1mcor","w2mag","w2sigm","w2flg",\
	"w2mcor","w3mag","w3sigm","w3flg","w3mcor","w4mag","w4sigm","w4flg","w4mcor","w1mag_1",\
	"w1sigm_1","w1flg_1","w2mag_1","w2sigm_1","w2flg_1","w3mag_1","w3sigm_1","w3flg_1",\
	"w4mag_1","w4sigm_1","w4flg_1","w1mag_2","w1sigm_2","w1flg_2","w2mag_2","w2sigm_2",\
	"w2flg_2","w3mag_2","w3sigm_2","w3flg_2","w4mag_2","w4sigm_2","w4flg_2","w1mag_3",\
	"w1sigm_3","w1flg_3","w2mag_3","w2sigm_3","w2flg_3","w3mag_3","w3sigm_3","w3flg_3",\
	"w4mag_3","w4sigm_3","w4flg_3","w1mag_4","w1sigm_4","w1flg_4","w2mag_4","w2sigm_4",\
	"w2flg_4","w3mag_4","w3sigm_4","w3flg_4","w4mag_4","w4sigm_4","w4flg_4","w1mag_5",\
	"w1sigm_5","w1flg_5","w2mag_5","w2sigm_5","w2flg_5","w3mag_5","w3sigm_5","w3flg_5",\
	"w4mag_5","w4sigm_5","w4flg_5","w1mag_6","w1sigm_6","w1flg_6","w2mag_6","w2sigm_6",\
	"w2flg_6","w3mag_6","w3sigm_6","w3flg_6","w4mag_6","w4sigm_6","w4flg_6","w1mag_7",\
	"w1sigm_7","w1flg_7","w2mag_7","w2sigm_7","w2flg_7","w3mag_7","w3sigm_7","w3flg_7",\
	"w4mag_7","w4sigm_7","w4flg_7","w1mag_8","w1sigm_8","w1flg_8","w2mag_8","w2sigm_8",\
	"w2flg_8","w3mag_8","w3sigm_8","w3flg_8","w4mag_8","w4sigm_8","w4flg_8","w1magp",\
	"w1sigp1","w1sigp2","w1k","w1ndf","w1mlq","w1mjdmin","w1mjdmax","w1mjdmean",\
	"w2magp","w2sigp1","w2sigp2","w2k","w2ndf","w2mlq","w2mjdmin","w2mjdmax","w2mjdmean",\
	"w3magp","w3sigp1","w3sigp2","w3k","w3ndf","w3mlq","w3mjdmin","w3mjdmax","w3mjdmean",\
	"w4magp","w4sigp1","w4sigp2","w4k","w4ndf","w4mlq","w4mjdmin","w4mjdmax","w4mjdmean",\
	"rho12","rho23","rho34","q12","q23","q34","xscprox","w1rsemi","w1ba","w1pa","w1gmag",\
	"w1gerr","w1gflg","w2rsemi","w2ba","w2pa","w2gmag","w2gerr","w2gflg","w3rsemi","w3ba",\
	"w3pa","w3gmag","w3gerr","w3gflg","w4rsemi","w4ba","w4pa","w4gmag","w4gerr","w4gflg",\
	"tmass_key","r_2mass","pa_2mass","n_2mass","j_m_2mass","j_msig_2mass","h_m_2mass",\
	"h_msig_2mass","k_m_2mass","k_msig_2mass","htm20"]

#Join longcols into a list so that you can dynamically feed newcols into the URL with the correct format.
newcols = ','.join(longcols) 

#Form the API POST object. Refer to these docs to select the desired parameters: https://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html 
# spatial -- spatial search area. Enter "Upload" to refer to a list of sources in your targets.tbl. 
# catalog -- name of the catalog to search (currently configured: allwise_p3as_psd)
# outfmt -- output format (currently configured for outfmt=6, XML output)
# selcols -- This is currently configured to receive newcols object, which is the correctly formatted longcols or "long form" option)

# Modify these parameters to change the API form POST:
spatial="Upload"
filename="@targets.tbl" #you must have a file called 'targets.tbl' saved in the same directory as this Python script, and it must be a valid IPAC table. You can validate your table here: https://irsa.ipac.caltech.edu/applications/TblValidator/
catalog="allwise_p3as_psd"
outfmt='6'
uradius=''

# Form the data payload:
data = {'spatial':spatial,'catalog':catalog,'uradius':uradius,'outfmt':outfmt,'selcols':newcols}

# Form the POST URL:
url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'

#POST the form data to the API using that URL, local targets.tbl file, and the payload, and store the response as r.
f = open('targets.tbl', 'rb')
r = requests.post(url, files={'filename':f,'Content-Type':'multipart/form-data'}, data=data)

#Write that data to a retrieved_data.text file
with open ('retrieved_data.xml', 'w') as output_file:
	output_file.write(r.text)

#Parse the XML to retrieve the output HTML URL for the object
tree = ET.parse('retrieved_data.xml')
root = tree.getroot()

# Grab the URL we need to view the table of this search:
for output in root.findall('output'):
	result_URL = output.find('outHTML').text

# Determine the controller of the user's Chrome browser
controller = webbrowser.get('chrome')

# Tell that controller to open a new window to the result data URL
controller.open_new(result_URL)
