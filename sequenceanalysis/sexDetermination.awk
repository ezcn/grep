BEGIN {
    xReads = 0
    yReads = 0
    autReads = 0

    xSites = 0
    ySites = 0
    autSites = 0
}
{
    chr = $1
    pos = $2
    cov = $3
    if(chr == "chrX") {
        xReads += cov
        xSites += 1
    }
    else if(chr == "chrY") {
        yReads += cov
        ySites += 1
    }
    else {
        autReads += cov
        autSites += 1
    }
}
END {
    OFS="\t"
    print("xCoverage", xSites > 0 ? xReads / xSites : 0)
    print("yCoverage", ySites > 0 ? yReads / ySites : 0)
    print("autCoverage", autSites > 0 ? autReads / autSites : 0)
}
