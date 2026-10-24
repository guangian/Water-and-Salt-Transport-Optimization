export class ApiError extends Error {
  code: string;
  status: number;
  constructor(message: string, code = 'unknown', status = 0) {
    super(message);
    this.name = 'ApiError';
    this.code = code;
    this.status = status;
  }
}

export function toApiError(e: unknown): ApiError {
  if (e instanceof ApiError) return e;
  if (e instanceof Error) return new ApiError(e.message || '请求失败', 'network', 0);
  return new ApiError(String(e), 'unknown', 0);
}
