import os


class ConfigParser:
    """
    A parser for molecular simulation configuration files.
    """

    def __init__(self, filepath):
        """
        Initialize the ConfigParser with the path to the configuration file.

        :param filepath: Path to the configuration file.
        """
        self.filepath = filepath
        self.parameters = {}
        self._set_default_values()
        self._parse_file()

    def _set_default_values(self):
        """
        Initialize default values for parameters.
        """
        self.parameters = {
            "MoleculeName": None,  # Required
            "ResidueName": None,  # Required
            "Charge": 0,  # Default charge
            "AdjustSymmetry": False,  # Default symmetry adjustment
            "AtomType": "gaff2",  # Default atom type
            "BoxType": "dodecahedron",  # Default box type
            "BoxEdgeDistance": 1.2,  # Default box edge distance in nm
            "BaseScratchDir": "/scratch/param/tmp",  # Default scratch directory
            "Cation": "NA",  # Default cation
            "Anion": "CL",  # Default anion
            "Multiplicity": 1,  # Default spin multiplicity
            "IonConcentration": 0.15,  # Default ion concentration for simulation
            "Protocol": "bcc",
            "CalculateGroupAverage": False,
            "IgnoreSymmetry": False,
        }

    def _parse_file(self):
        """
        Parse the configuration file and populate the parameters dictionary.
        """
        if not os.path.isfile(self.filepath):
            raise FileNotFoundError(f"Configuration file '{self.filepath}' not found.")

        with open(self.filepath, "r") as file:
            for line_num, line in enumerate(file, start=1):
                stripped_line = line.strip()

                # Skip empty lines and full-line comments
                if not stripped_line or stripped_line.startswith("#"):
                    continue

                # Remove inline comments
                if "#" in stripped_line:
                    stripped_line = stripped_line.split("#", 1)[0].strip()

                # Split the line into key and value
                if ":" in stripped_line:
                    key, value = stripped_line.split(":", 1)
                    key = key.strip()
                    value = value.strip()

                    # Assign the parameter with appropriate type conversion
                    self._assign_parameter(key, value, line_num)
                else:
                    print(
                        f"Warning: Line {line_num} in '{self.filepath}' is not a valid key-value pair and will be ignored."
                    )

        # After parsing, validate required parameters
        self._validate_required_parameters()

    def _assign_parameter(self, key, value, line_num):
        """
        Assign the value to the corresponding key in parameters with type conversion.

        :param key: The parameter name.
        :param value: The parameter value as a string.
        :param line_num: The line number in the configuration file (for error reporting).
        """
        # Mapping of keys to their exact parameter names
        key_map = {
            "moleculename": "MoleculeName",
            "residuename": "ResidueName",
            "charge": "Charge",
            "adjustsymmetry": "AdjustSymmetry",
            "atomtype": "AtomType",
            "boxtype": "BoxType",
            "boxedgedistance": "BoxEdgeDistance",
            "basescratchdir": "BaseScratchDir",
            "cation": "Cation",
            "anion": "Anion",
            "multiplicity": "Multiplicity",
            "ionconcentration": "IonConcentration",
            "protocol": "Protocol",
            "calculategroupaverage": "CalculateGroupAverage",
            "ignoresymmetry": "IgnoreSymmetry",
        }

        key_lower = key.lower()

        if key_lower in key_map:
            actual_key = key_map[key_lower]
            try:
                if actual_key == "Charge":
                    self.parameters[actual_key] = int(value)
                elif actual_key == "AdjustSymmetry":
                    self.parameters[actual_key] = self._str_to_bool(value)
                elif actual_key == "IgnoreSymmetry":
                    self.parameters[actual_key] = self._str_to_bool(value)
                elif actual_key == "CalculateGroupAverage":
                    self.parameters[actual_key] = self._str_to_bool(value)
                elif actual_key == "BoxEdgeDistance":
                    self.parameters[actual_key] = float(value)
                elif actual_key == "Multiplicity":
                    self.parameters[actual_key] = int(value)
                elif actual_key == "IonConcentration":
                    self.parameters[actual_key] = float(value)
                else:
                    self.parameters[actual_key] = value
            except ValueError as ve:
                raise ValueError(
                    f"Invalid value for '{key}' on line {line_num}: '{value}'. Expected type {self._expected_type(actual_key).__name__}."
                ) from ve
        else:
            print(
                f"Warning: Unknown parameter '{key}' on line {line_num}. It will be ignored."
            )

    def _str_to_bool(self, value):
        """
        Convert a string to a boolean.

        :param value: The string to convert.
        :return: Boolean value.
        """
        value_lower = value.lower()
        if value_lower in ["true", "yes", "1", "TRUE"]:
            return True
        elif value_lower in ["false", "no", "0", "FALSE"]:
            return False
        else:
            raise ValueError(f"Cannot convert '{value}' to boolean.")

    def _expected_type(self, key):
        """
        Return the expected type for a given parameter key.

        :param key: The parameter name.
        :return: The expected Python type.
        """
        type_map = {
            "Charge": int,
            "AdjustSymmetry": bool,
            "BoxEdgeDistance": float,
            "Multiplicity": int,
            # Add more mappings if needed
        }
        return type_map.get(key, str)

    def _validate_required_parameters(self):
        """
        Validate that all required parameters are present.
        """
        required = ["MoleculeName", "ResidueName"]
        missing = [param for param in required if not self.parameters.get(param)]
        if missing:
            raise ValueError(f"Missing required parameters: {', '.join(missing)}")

    def get_parameter(self, key):
        """
        Get the value of a parameter by key.

        :param key: The parameter name.
        :return: The parameter value.
        """
        return self.parameters.get(key)

    def __getattr__(self, name):
        """
        Allow attribute-style access to parameters.

        :param name: The parameter name.
        :return: The parameter value.
        """
        if name in self.parameters:
            return self.parameters[name]
        raise AttributeError(f"'ConfigParser' object has no attribute '{name}'")

    def __repr__(self):
        """
        String representation of the ConfigParser instance.

        :return: String representation.
        """
        return f"ConfigParser({self.parameters})"
