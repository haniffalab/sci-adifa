class InternalServerError(Exception):
    pass


class SchemaValidationError(Exception):
    pass


class DatabaseOperationError(Exception):
    pass


class InvalidDatasetIdError(Exception):
    pass


class InvalidModalityError(Exception):
    pass


class DatasetNotExistsError(Exception):
    pass


class UnauthorizedError(Exception):
    pass


errors = {
    "InternalServerError": {"message": "Something went wrong", "status": 500},
    "SchemaValidationError": {
        "message": "Request is missing required fields",
        "status": 400,
    },
    "DatabaseOperationError": {"message": "DatabaseOperationError", "status": 400},
    "InvalidDatasetIdError": {
        "message": "Invalid Dataset ID passed in request",
        "status": 400,
    },
    "DatasetNotExistsError": {
        "message": "Dataset with given id doesn't exist",
        "status": 400,
    },
    "SchemaValidationError": {
        "message": "Request is missing required fields",
        "status": 400,
    },
    "DatabaseOperationError": {"message": "DatabaseOperationError", "status": 400},
    "InvalidDatasetIdError": {
        "message": "Invalid Dataset ID passed in request",
        "status": 400,
    },
    "InvalidModalityError": {
        "message": "Invalid modality passed in request",
        "status": 400,
    },
    "DatasetNotExistsError": {
        "message": "Dataset with given id doesn't exist",
        "status": 400,
    },
    "UnauthorizedError": {"message": "Invalid username or password", "status": 401},
}
