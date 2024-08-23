'''Parsekey wrangling module. Adapted from Emily Bentley's nanodrop wrangling utilities'''

class ParseKey: #Taken from EB NanoDrop wrangling file and modified by Jenna T. to accommodate lists
    """ParseKeys for use in _make_columns_by_parse_key and analyze_sample_names."""
    
    def __init__(self, *args, separator=None):
        """
        Initialize a ParseKey to interpret the "Sample ID" column of a DataFrame.
            
        *args : any number of two-tuples containing (column name, type) pairs
            these must be provided in the same order they appear in the "Sample ID" column
            all column names must be unique
        separator : str, default None
            distinguishes separate elements in the "Sample ID" column
            must be specified by kwarg; default is only allowed for a single arg arg
        
        For example:
            parse_rna_peptide = ParseKey(
                ("peptide", str),
                ("concentration", float),
                ("ratio", float),
                separator="_",)
            would be use to interpret a "Sample ID" column whose contents take the form
            peptide_concentration_ratio
        
        handle_input.request_parsekey_specifications() is a command-line utility that will assist in building a ParseKey.
        """
        
        if separator == None:
            try:
                assert len(args) == 1
            except AssertionError:
                raise RuntimeError("A separator must be specified by kwarg when more than one arg is given.")
        
        for arg in args:
            try:
                assert type(arg) == tuple
                assert len(arg) == 2
                assert type(arg[0]) == str or type(arg[0][0]) == str #second part is to handle lists: ([name,size],[str,float])
                assert type(arg[1]) == type or type(arg[1][0]) == type #second part is to handle lists: ([name,size],[str,float])
            except AssertionError:
                raise RuntimeError(f"\n{arg} is not a valid input."
                                   "Expected a tuple of length 2 where item 0 is a str and item 1 is a type.")
        
        self.separator = separator
        self.parse_key = tuple([*args])
        
class ParseName:
    def __init__(self,string_to_split, parsekey):
        self.string = string_to_split
        self.parsekey = parsekey
        self.parsekey_dict = {}
    
    def parse_name(self):
        parsekey = self.parsekey
        fn = self.string
        if isinstance(parsekey.separator,list) == True:
            main_sep = parsekey.separator[0]
            sub_sep = parsekey.separator[1]
        else:
            main_sep = parsekey.separator
        
        index_main_sep = 0
        filename_parse_list = fn.split(main_sep)
        
        try:
            assert len(filename_parse_list) == len(parsekey.parse_key)
        except AssertionError:
            raise RuntimeError("Parsekey does not have the same number of entries as the filepath.")
            
        for key,datatype in parsekey.parse_key:
#             print(key, datatype)
            if isinstance(datatype,list) == True:
                second_sep = filename_parse_list[index_main_sep].split(sub_sep)
                index_subsep = 0
                if isinstance(key,list) == True:
                    try:
                        assert len(key) == len(datatype)
                    except AssertionError:
                        raise RuntimeError(f"\n{key} and {datatype} is not a valid input."
                                   "Expected lists of equal length.")
                    for item in key:
                        if datatype[index_subsep] == float or datatype[index_subsep] == int:
                            val = self.findunit(second_sep[index_subsep])[1]
                            unit = self.findunit(second_sep[index_subsep])[0]
                            if unit == "None":
                            	cat_name = item
                            else:
                            	cat_name = f'{item} ({unit})'
                        else:
#                             print(item,datatype[index_subsep])
                            val = second_sep[index_subsep]
                            cat_name = item
                        self.parsekey_dict[cat_name] = val
                        index_subsep +=1
                else:
                    cat_name = key
                    conc_unit = self.findunit(second_sep[1])[0]
                    conc_value = self.findunit(second_sep[1])[1]
                    cat_conc = f'[{cat_name}] ({conc_unit})'
                    self.parsekey_dict[cat_name] = second_sep[0]
                    self.parsekey_dict[cat_conc] = conc_value
            else:
                name = filename_parse_list[index_main_sep]
                if ".txt" in name:
                	name = name.replace(".txt","")
                if '.tif' in name:
                    name = name.replace(".tif","")
                if "FRAP" in name:
                    name = name.replace("FRAP","")
                    if "T" in name or "S" in name or "U" in name:
                        name = name[1:]
                    name = int(name)
                if datatype == float or datatype == int:
                	unit, val = self.findunit(name)
                	key = f'{key} ({unit})'
                	name = val
                
                self.parsekey_dict[key] = name
            
            index_main_sep += 1

    def findunit(self,test_string):
        '''Takes a string and finds the unit indicated. Default unit = mg/mL, returns tuple of "unit", value (float)'''
        if ".txt" in test_string:
            test_string = test_string.replace(".txt","")
        elif ".csv" in test_string:
        	test_string = test_string.replace(".csv","")
        unit = "None"
        unit_val = False
        i = -1
        while unit_val == False:
            if test_string[i].isdigit() == True or test_string[i] == ".":
                unit_val = True
                if i < -1:
                    unit = str(test_string[i+1:])
                    number = float(test_string[:i+1])
                    if len(unit) == 4 and unit.lower()[1] == "g" and unit.lower()[-1]=="l":
                        old_str = unit.lower()
                        unit = old_str[0:2] + "/" + old_str[2] + "L"
                    elif unit=="p": 
                    	unit = "percent"
                else:
                    number = float(test_string)
            else:
                i = i - 1  
        return (unit,number)
    
    def return_parsekey_dict(self):
        return self.parsekey_dict
