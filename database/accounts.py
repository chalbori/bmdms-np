class Database:
    # The public BMDMS-NP server is no longer reachable from outside the BMDRC network.
    # Set these values to a MySQL/MariaDB database that you host yourself.
    host = 'localhost'
    database = 'bmdms-np'
    user_id = 'your_user_id'
    user_password = 'your_password'

    # @property
    # def host(self):
    #     return self._host
    #
    # @property
    # def database(self):
    #     return self._database
    #
    # @property
    # def user_id(self):
    #     return self._user_id
    #
    # @property
    # def user_password(self):
    #     return self._user_password
