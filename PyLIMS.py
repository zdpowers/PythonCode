import numpy as np
import pandas as pd
import itertools

class Rack:
  
  #positions = 24

  def __init__(self, id):
    self.id = id
    self.contents = []

  def add_sample(self, sample):
    if len(self.contents) < 24:
      self.contents.append(sample)
    else:
      print("Rack is full. Add sample to new rack.")
  
  def add_sample_group(self, samplelist):
    qty = len(samplelist)
    openings = 24 - len(self.contents)
    if qty > openings:
      added = samplelist[:openings]
      remaining = samplelist[openings:]
      self.contents.extend(added)
      print("Samples Not Added: " + str(remaining))
    elif qty <= openings:
      self.contents.extend(samplelist)

  def remove_sample(self, sampleid):
    position = self.contents.index(sampleid)
    self.contents[position] = ""
    print(f"Sample {sampleid} removed at position {position}")

  def add_at_index(self, sample, position):
    current_position_value = self.contents[position]
    if  current_position_value == "":
      self.contents[position] = sample
    else:
      print(f"Position {position} already contains {current_position_value}")

  def replace_at_index(self, sample, position):
    current_position_value = self.contents[position]
    self.contents[position] = sample
    print(f"Sample {sample} replaced {current_position_value} at position {position}")

  def get_empty_positions(self):
    indices = [index for index, element in enumerate(self.contents) if element == '']
    print(indices)

  def get_sample_list(self):
    samplelist = []
    for i in self.contents:
      if i != '':
        samplelist.append(i)
    print(samplelist)

  def get_sample_count(self):
    samplelist = []
    for i in self.contents:
      if i != '':
        samplelist.append(i)
    print(len(samplelist))

  def rack_array(self):
      return(np.array(self.contents))
	
	
class Plate:

  positions = 96
  columns = list(map(lambda x: chr(ord('@')+ x), range(1,13)))
  rows = list(range(1,13))

  def __init__(self, id):
    self.id = id
    self.racks = []
    self.cols = []
    self.contents = []
    self.columns = {}

  def add_rack(self, rack):
    qty = len(rack) # Quantity of samples on rack
    rck = rack
    column_count = len(self.columns) # Number of columns in the column dicitonary
    # If the rack has less then 24 samples, fill out the rest of the rack with blanks
    if qty < 24: 
      diff = 24 - qty
      for x in range(0, diff):
        rck.append('')
        #rck.append(np.nan)

    openings = 96 - len(self.contents) # Number of wells left on the plate
    quantity = len(rck) # number of samples in the rack after it has been filled out with blanks

    if quantity <= openings:
      self.racks.append(rck) # Add the rack to the list of racks
      self.contents.extend(rck)  # Add samples to list of the plates samples
      
      self.cols.append(rck[0:8]) # Slice first set of 8 samples on rack and add them to the column list
      self.columns[len(self.cols)] = pd.Series(rck[0:8], index=list(map(lambda x: chr(ord('@')+ x), range(1,9))))

      self.cols.append(rck[8:16]) # Slice second set of 8 of samples on rack and add them to the column list
      self.columns[len(self.cols)] = pd.Series(rck[8:16], index=list(map(lambda x: chr(ord('@')+ x), range(1,9))))

      self.cols.append(rck[16:24]) # Slice third set of 8 of samples on rack and add them to the column list
      self.columns[len(self.cols)] = pd.Series(rck[16:24], index=list(map(lambda x: chr(ord('@')+ x), range(1,9))))

  def add_racks(self, racks):
    # TODO: Function to add list of racks
    pass


  def plate_df(self):
    column_count = len(self.columns) # Number of columns in the column dicitonary
    remaining_cols = 12 - column_count # Number of columns in the column dicitonary
    blanks = [] # List of blanks
    for i in range(0,8):
      blanks.append("")
      #blanks.append(np.nan)
    # first fill out the rest of the plate with blanks
    if column_count < 12:
      for i in range(0, remaining_cols):
        self.cols.append(blanks)
        self.columns[len(self.cols)] = pd.Series(blanks, index=list(map(lambda x: chr(ord('@')+ x), range(1,9))))
    
    # Create df
    df = pd.DataFrame(self.columns)
    return(df)

class PCRPlate:

  def __init__(self, id):
    self.id = id
    self.plates = []
    self.contents = []
    self.quadrants = {}

  def add_plate_df(self, plate):
    if len(self.plates) < 4:
      self.plates.append(plate)
      self.quadrants[len(self.plates)] = plate
  
  def fill_remaining(self):
    num_plates = len(self.plates)
    quads_remaining = 4 - num_plates
    blank_dict = {}
    blank_list = []
    blanks = [] # List of blanks
    for i in range(0,8):
      blanks.append("")
    
    for i in range(0, 12):
      blank_list.append(blanks)  
      blank_dict[len(blank_list)] = pd.Series(blanks, index=list(map(lambda x: chr(ord('@')+ x), range(1,9))))

    
    empty_quad = pd.DataFrame(blank_dict)
    if num_plates < 4:
      for i in range(0, quads_remaining):
        self.plates.append(empty_quad)
        self.quadrants[len(self.plates)] = empty_quad

  # The following 2 work, ignore the fill remaining function
  def make_blank(self):
    blankplate = []
    for i in range(0,12):
      blank_col = ["","","","","","","",""]
      blankplate.append(blank_col)
    return np.array(blankplate)

  def finishplate(self):
    num_plates = len(self.plates)
    quads_remaining = 4 - num_plates
    if num_plates < 4:
      for i in range(0, quads_remaining):
          self.plates.append(self.make_blank())

  def make_plate(self):
    cols1 = np.split(self.plates[0], 12)
    cols2 = np.split(self.plates[1], 12)
    cols3 = np.split(self.plates[2], 12)
    cols4 = np.split(self.plates[3], 12)
    q1q3 = [x for x in itertools.chain.from_iterable(itertools.zip_longest(cols1, cols3))]
    q2q4 = [x for x in itertools.chain.from_iterable(itertools.zip_longest(cols2, cols4))]
    p1p3 = np.vstack(q1q3).T
    p2p4 = np.vstack(q2q4).T
    q13 = np.split(p1p3, 8)
    q24 = np.split(p2p4, 8)
    pcr1234 = [x for x in itertools.chain.from_iterable(itertools.zip_longest(q13, q24))]
    pcr = np.vstack(pcr1234)
    return pcr
  
