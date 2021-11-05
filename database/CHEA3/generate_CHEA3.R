# Script to download CHEA3 data

# Download the latest version of CHEA3 ENCODE data
file = "./ENCODE_ChIP-seq.gmt"
url <- "https://maayanlab.cloud/chea3/assets/tflibs/ENCODE_ChIP-seq.gmt"

curl::curl_download(url, file)

# Download the latest version of CHEA3 REMAP data
file = "./ReMap_ChIP-seq.gmt"
url <- "https://maayanlab.cloud/chea3/assets/tflibs/ReMap_ChIP-seq.gmt"

curl::curl_download(url, file)

# Download the latest version of CHEA3 Literature data
file = "./Literature_ChIP-seq.gmt"
url <- "https://maayanlab.cloud/chea3/assets/tflibs/Literature_ChIP-seq.gmt"

curl::curl_download(url, file)