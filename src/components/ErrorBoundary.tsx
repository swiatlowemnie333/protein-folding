import React, { Component, ErrorInfo, ReactNode } from 'react';
import { AlertTriangle, RotateCcw } from 'lucide-react';

interface Props {
  children: ReactNode;
  fallbackTitle?: string;
  onReset?: () => void;
}

interface State {
  hasError: boolean;
  error: Error | null;
}

export class ErrorBoundary extends Component<Props, State> {
  constructor(props: Props) {
    super(props);
    this.state = {
      hasError: false,
      error: null
    };
  }

  public static getDerivedStateFromError(error: Error): State {
    return { hasError: true, error };
  }

  public componentDidCatch(error: Error, errorInfo: ErrorInfo) {
    console.error('ErrorBoundary caught an error:', error, errorInfo);
  }

  private handleReset = () => {
    this.setState({ hasError: false, error: null });
    if (this.props.onReset) {
      this.props.onReset();
    }
  };

  public render(): ReactNode {
    if (this.state.hasError) {
      return (
        <div className="flex-1 w-full h-full min-h-[300px] flex flex-col items-center justify-center bg-black/95 text-white p-6 font-mono select-none">
          <div className="max-w-md w-full bg-zinc-950 border border-red-800/80 rounded-xl p-5 shadow-2xl flex flex-col items-center text-center space-y-4">
            <div className="w-12 h-12 rounded-full bg-red-950/80 border border-red-500/50 flex items-center justify-center text-red-400">
              <AlertTriangle className="w-6 h-6 text-red-400" />
            </div>
            
            <div className="space-y-1">
              <h3 className="text-sm font-bold text-white uppercase tracking-wider">
                {this.props.fallbackTitle || 'Wykryto problem z renderowaniem 3D'}
              </h3>
              <p className="text-xs text-zinc-400">
                Nastąpiło tymczasowe rozsynchronizowanie buforów geometrii. Kliknij poniżej, aby zresetować stan widoku.
              </p>
            </div>

            {this.state.error && (
              <div className="w-full bg-black/80 border border-zinc-800 rounded p-2 text-[10px] text-red-300 text-left font-mono overflow-x-auto max-h-24">
                {this.state.error.message}
              </div>
            )}

            <button
              onClick={this.handleReset}
              className="px-4 py-2 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 text-white text-xs font-bold rounded-lg flex items-center gap-2 cursor-pointer shadow-lg transition"
            >
              <RotateCcw className="w-3.5 h-3.5" />
              <span>Zresetuj widok 3D</span>
            </button>
          </div>
        </div>
      );
    }

    return this.props.children;
  }
}
