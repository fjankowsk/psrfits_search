# 				RESCALE TIME
# 	Resize data array extracted from a FITS file to a new number of rows.
# 	Used to decrease the time lapse of a block for a involve the rfifind processing.


# MODULES

import numpy as np
import astropy.io.fits as fi
import sys
import argparse as arg
import os


# ARGUMENTS LIST

parser = arg.ArgumentParser(
    description="Resize the FITS file to decrease the time lapse of a row."
)

parser.add_argument("fileName", type=str, help="Name of the FITS file to adjust.")
parser.add_argument(
    "-n",
    dest="n",
    type=float,
    required=True,
    help="Dividing factor on samples per block. Equivalent as multiplying factor for the number of row.",
)
parser.add_argument(
    "-o",
    dest="newFileName",
    type=str,
    required=True,
    help="Name of the new FITS file to write.",
)

args = parser.parse_args()


# CHECKING INPUT PARAMETERS

if os.path.isfile(args.fileName):  # Checking file existence
    print("\nExtraction of data from {:s}.\n".format(args.fileName))
else:
    print("\n{:s} is not a file.\n".format(args.fileName))
    sys.exit()

print("Resized arrays writed in {:s}.\n".format(args.newFileName))


# DATA EXTRACTION OF THE PREVIOUS FITS

hdul = fi.open(args.fileName, mode="readonly")

headObs = hdul[0].header
head = hdul[1].header
data = hdul[1].data

chan = head["NCHAN"]  # Number of frequency channel
samples = head["NSBLK"]  # Number of samples per block
blocks = head["NAXIS2"]  # Number of blocks
pol = head["NPOL"]  # Number of polarizations
bin = head["NBIN"]  # Number of bins ( normally 1 )
bits = head["NBITS"]  # Number of bits by datum
tsample = head["TBIN"]  # Time lapse of a sample
lst = headObs["STT_LST"]

print("Old samples per block: {0}".format(samples))
print("Old blocks: {0}".format(blocks))
print("Old block size: {0}".format(bin * chan * pol * samples))
print("Old time per block: {0}".format(tsample * samples))

# COMPUTING OF THE NEW INFORMATION

newSamples = int(samples / args.n)  # Computing of the new number of samples per block
newBlocks = int(blocks * args.n)  # Computing of the new number of blocks
newSize = int(
    bin * chan * pol * newSamples
)  # Computing of the new size of the data array
newTblock = float(tsample * newSamples)  # Computing of the new block time lapse
newLst = round(lst - tsample * samples / 2.0 + newTblock / 2.0, 0)

print("New samples per block: {0}".format(newSamples))
print("New blocks: {0}".format(newBlocks))
print("New block size: {0}".format(newSize))
print("New time per block: {0}".format(newTblock))

# WRITING MODIFIED LINES IN HEADERS

head["NAXIS1"] = newSize  # new data array size (block size)
head["NAXIS2"] = newBlocks  # new number of rows in the data fits
head["TFORM17"] = str(newSize) + "E"  # new amplitude data array byte size
head["NSBLK"] = newSamples  # number of samples per block in the new fits file
head["TDIM17"] = (
    "(" + str(bin) + "," + str(chan) + "," + str(pol) + "," + str(newSamples) + ")"
)  # new amplitude data array dimension
headObs["STT_LST"] = newLst  # new LST of the first block center

if (
    newSize < bits
):  # Checking if there are a number of samples greater than the number of bits
    print("The new number of samples per block is less than {:d}.".format(bits))
    print(
        "The minimal samples number is {:d}. Thus, n must be less or equal than {:d}.\n".format(
            bits, samples / bits
        )
    )
    sys.exit()


# RESIZING ARRAYS
cols = hdul[1].columns
print(cols)
print(cols.names)

# tsubint (time per block)
newArray = np.resize(
    data["tsubint"] / args.n, (newBlocks,)
)  # Computing of the new values and resizing of the data array
print(newArray.shape)
data["tsubint"] = newArray

# offs_sub
newArray = np.arange(
    newTblock / 2.0, tsample * samples * blocks, newTblock
)  # Creation of the new block offset 1D array
print(newArray.shape)
data["offs_sub"] = newArray

# lst_sub
newArray = np.around(
    np.linspace(newLst, newLst + newBlocks * newTblock, newBlocks), 0
)  # Creation of the new LST time 1D array
print(newArray.shape)
data["lst_sub"] = newArray

for f in range(3, 12):  # Loop on other 1D subint arrays
    field = cols.names[f]
    newArray = np.resize(data[field], (newBlocks,))  # Resizing of the data array
    print(field, newArray.shape)
    data[field] = newArray

for f in range(12, 14):  # Loop on 2D weight arrays
    field = cols.names[f]
    newArray = np.resize(data[field], (newBlocks, chan))  # Resizing of the data array
    print(field, newArray.shape)
    data[field] = newArray

for f in range(14, 16):  # Loop on 2D weight arrays
    field = cols.names[f]
    newArray = np.resize(
        data["field"], (newBlocks, chan * pol)
    )  # Resizing of the data array
    print(field, newArray.shape)
    data[field] = newArray

# data
newFormat = fi.column._ColumnFormat(
    str(newSize) + "E"
)  # Definition of the new data array format
newDim = (
    "(" + str(bin) + "," + str(chan) + "," + str(pol) + "," + str(newSamples) + ")"
)  # Definition of the new data array definition
newArray = np.reshape(
    data["data"], (newBlocks, newSamples, pol, chan, bin)
)  # Resizing of the data array
print(newArray.shape)
data["data"] = newArray
cols.change_attrib("data", "format", newFormat)
cols.change_attrib("data", "dim", newDim)

# DEFINITION OF THE NEW FITS

print(head)

hdul.writeto(args.newFileName)
