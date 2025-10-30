# users

Methods to interact with the user info on the API.

### get_costs(start_year: int, start_month: int, end_year: int | None = None, end_month: int | None = None) → List[CostDetail]

Get the user costs.

This function gets a dict with the user costs.

* **Returns:**
  Dict with the user costs.

### get_fees() → UserFees

Get information about the user’s fees.

### get_info() → User

Get the user information.

This funtion gets the user information, including the user’s name, email,
username, plan, programs, and total available credits.

* **Returns:**
  Dict with the user information.

### get_quotas() → Dict[str, Quota]

Get the user quotas.

This function gets a dict with the user quotas.

* **Returns:**
  Dict with the user quotas.

### get_task_orchestration_fee_warning()

::docsbannersmall
::
