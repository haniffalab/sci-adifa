from werkzeug.exceptions import HTTPException


class ApiException(HTTPException):
    code = 500
    description = 'Internal Server Error'
    def __init__(self, e, description=None, response=None):
        if description is not None:
            self.description = description        
        self.errors = e


class InternalServerError(ApiException):
    pass


class SqlAlchemyError(ApiException):
    code = 500
    description = 'Database Error'


class SchemaValidationError(ApiException):
    code = 400
    description = 'Request is missing required fields'


class EmailAlreadyExistsError(ApiException):
    code = 400
    description = 'User with given email address already exists'


class UnauthorizedError(ApiException):
    code = 401
    description = 'Invalid username or password'


class EmailDoesnotExistsError(ApiException):
    code = 400
    description = 'Couldn\'t find the user with given email address'


class BadTokenError(ApiException):
    code = 403
    description = 'Invalid token'
